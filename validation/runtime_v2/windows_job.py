"""Windows-only owned process trees. No taskkill / image-name termination.

Uses the audited frozen structure definitions, but not its fixed four-core
runner. Child handle inheritance is restricted to stdin/stdout, and the
supervisor owns the Job handle throughout. No restricted startup affinity:
Gaussian 16W A.03's PGI runtime is known to fault in that configuration.
"""
from __future__ import annotations

import ctypes as ct
from ctypes import wintypes as wt
import os
import json
from pathlib import Path
import subprocess
import threading
import time

from fast_checkpoint import ActiveClock
from windows_pause import JobPause

SUPPORTS_RAM_PAUSE = True


def api():
    if os.name != 'nt' or ct.sizeof(ct.c_void_p) != 8:
        raise RuntimeError('Production execution requires 64-bit Windows Python')
    import validation_process as v
    k = v.kernel()
    signatures = {
        'CreateJobObjectW': ([ct.c_void_p, wt.LPCWSTR], wt.HANDLE),
        'SetInformationJobObject': ([wt.HANDLE, ct.c_int, ct.c_void_p, wt.DWORD], wt.BOOL),
        'QueryInformationJobObject': ([wt.HANDLE, ct.c_int, ct.c_void_p, wt.DWORD, ct.c_void_p], wt.BOOL),
        'AssignProcessToJobObject': ([wt.HANDLE, wt.HANDLE], wt.BOOL),
        'TerminateJobObject': ([wt.HANDLE, wt.UINT], wt.BOOL),
        'TerminateProcess': ([wt.HANDLE, wt.UINT], wt.BOOL),
        'GetExitCodeProcess': ([wt.HANDLE, ct.POINTER(wt.DWORD)], wt.BOOL),
        'CloseHandle': ([wt.HANDLE], wt.BOOL),
        'ResumeThread': ([wt.HANDLE], wt.DWORD),
        'WaitForSingleObject': ([wt.HANDLE, wt.DWORD], wt.DWORD),
        'GetCurrentProcess': ([], wt.HANDLE),
        'GetProcessAffinityMask': ([wt.HANDLE, ct.POINTER(ct.c_size_t), ct.POINTER(ct.c_size_t)], wt.BOOL),
        'SetProcessAffinityMask': ([wt.HANDLE, ct.c_size_t], wt.BOOL),
        'GetActiveProcessorGroupCount': ([], wt.WORD),
        'InitializeProcThreadAttributeList': ([ct.c_void_p, wt.DWORD, wt.DWORD, ct.POINTER(ct.c_size_t)], wt.BOOL),
        'UpdateProcThreadAttribute': ([ct.c_void_p, wt.DWORD, ct.c_size_t, ct.c_void_p, ct.c_size_t, ct.c_void_p, ct.c_void_p], wt.BOOL),
        'DeleteProcThreadAttributeList': ([ct.c_void_p], None),
        'CreateProcessW': ([wt.LPCWSTR, wt.LPWSTR, ct.c_void_p, ct.c_void_p, wt.BOOL, wt.DWORD,
                            ct.c_void_p, wt.LPCWSTR, ct.c_void_p, ct.POINTER(v.PROCESSINFO)], wt.BOOL),
    }
    for name, (args, result) in signatures.items():
        fun = getattr(k, name)
        fun.argtypes, fun.restype = args, result
    return k, v


def host_info() -> dict:
    k, v = api()
    if k.GetActiveProcessorGroupCount() != 1:
        raise RuntimeError('Multiple processor groups require a reviewed topology policy; refusing to miscount')
    masks = v.physical_core_masks(64)
    class MEMORY(ct.Structure):
        _fields_ = [('length', wt.DWORD), ('load', wt.DWORD)] + [(x, ct.c_uint64) for x in
            ('total', 'available', 'page_total', 'page_available', 'virtual_total', 'virtual_available', 'extended')]
    memory = MEMORY()
    memory.length = ct.sizeof(memory)
    k.GlobalMemoryStatusEx.argtypes = [ct.POINTER(MEMORY)]
    v.require(k.GlobalMemoryStatusEx(ct.byref(memory)))
    return {'physical_cores': len(masks), 'logical_processors': os.cpu_count(),
            'physical_masks': masks, 'total_gib': memory.total / 2**30,
            'available_gib': memory.available / 2**30,
            'profile': '98x3-small-host-policy' if len(masks) < 16 else '16-plus-core-policy',
            'topology_source': 'Windows GetLogicalProcessorInformation; one record per physical core',
            'affinity_policy': 'full startup affinity; threads and Job CPU-rate cap, not exclusive core pinning'}


def full_startup_affinity() -> None:
    k, v = api()
    current, system = ct.c_size_t(), ct.c_size_t()
    process = k.GetCurrentProcess()
    v.require(k.GetProcessAffinityMask(process, ct.byref(current), ct.byref(system)))
    v.require(k.SetProcessAffinityMask(process, system.value))


def existing_gaussian(gaussian: Path) -> list[dict]:
    """Read-only admission check. Do not kill a foreign Gaussian/legacy runner."""
    import re
    class Process(ct.Structure):
        _fields_ = [('size', wt.DWORD), ('usage', wt.DWORD), ('pid', wt.DWORD), ('heap', ct.c_size_t),
            ('module', wt.DWORD), ('threads', wt.DWORD), ('parent', wt.DWORD), ('priority', wt.LONG),
            ('flags', wt.DWORD), ('name', wt.WCHAR * 260)]
    k, _ = api()
    k.CreateToolhelp32Snapshot.argtypes = [wt.DWORD, wt.DWORD]
    k.CreateToolhelp32Snapshot.restype = wt.HANDLE
    for name in ('Process32FirstW', 'Process32NextW'):
        getattr(k, name).argtypes = [wt.HANDLE, ct.POINTER(Process)]
        getattr(k, name).restype = wt.BOOL
    snapshot = k.CreateToolhelp32Snapshot(2, 0)
    if snapshot == ct.c_void_p(-1).value:
        raise ct.WinError(ct.get_last_error())
    records = []
    deadline = time.monotonic() + 10
    try:
        process = Process(); process.size = ct.sizeof(process)
        present = k.Process32FirstW(snapshot, ct.byref(process))
        while present:
            if time.monotonic() >= deadline:
                raise TimeoutError('Process admission snapshot exceeded its deadline')
            if re.fullmatch(r'(g16|l[0-9]+)\.exe', process.name, re.I):
                records.append({'ProcessId': process.pid, 'Name': process.name})
            present = k.Process32NextW(snapshot, ct.byref(process))
        if ct.get_last_error() != 18:  # ERROR_NO_MORE_FILES
            raise ct.WinError(ct.get_last_error())
    finally:
        k.CloseHandle(snapshot)
    return records


_LAUNCH_LOCK = threading.Lock()


def _drain_owned_tree(k, v, job, timeout=5.0):
    accounting = v.ACCOUNTING()
    v.require(k.QueryInformationJobObject(job, 1, ct.byref(accounting), ct.sizeof(accounting), None))
    if accounting.ActiveProcesses:
        v.require(k.TerminateJobObject(job, 125))
    deadline = time.monotonic() + timeout
    while accounting.ActiveProcesses:
        if time.monotonic() >= deadline:
            raise RuntimeError('Owned Job cleanup did not confirm zero processes before its deadline')
        time.sleep(.05)
        v.require(k.QueryInformationJobObject(job, 1, ct.byref(accounting), ct.sizeof(accounting), None))


def run_tree(argv: list, cwd: Path, cores: int, memory_gib: int, timeout: float,
             cancel=lambda: False, on_started=lambda x: None, on_tick=lambda x: None,
             hold=lambda: False, on_hold_error=lambda x: None,
             on_timeout=None) -> dict:
    if not 1 <= cores <= 14 or memory_gib < 1 or timeout <= 0:
        raise ValueError('Invalid resource envelope')
    k, v = api()
    import msvcrt
    class STARTUPINFOEX(ct.Structure):
        _fields_ = [('StartupInfo', v.STARTUPINFO), ('lpAttributeList', ct.c_void_p)]
    job = v.require(k.CreateJobObjectW(None, None))
    proc = v.PROCESSINFO()
    attributes = None
    initialized = False
    stdin = console = None
    begin, started = time.monotonic(), time.time()
    reason = None
    pause = None
    clock = ActiveClock(begin)
    timeout_notified = False
    timeout_count = 0
    next_timeout = timeout
    previous_wanted = False
    termination_started = None
    hold_failed = False
    pause_error = None
    pause_failures = []
    last_pause_scan = 0.0
    try:
        pause = JobPause(k, job)
        stdin = open(os.devnull, 'rb')
        console = (cwd / 'launcher.log').open('ab')
        limits = v.EXTENDEDLIMIT()
        limits.BasicLimitInformation.LimitFlags = 0x200 | 0x2000  # JOB_MEMORY | KILL_ON_JOB_CLOSE
        limits.JobMemoryLimit = memory_gib * 2**30
        v.require(k.SetInformationJobObject(job, 9, ct.byref(limits), ct.sizeof(limits)))
        rate_value = max(1, min(10000, int(10000 * cores / os.cpu_count())))
        rate = (wt.DWORD * 2)(0x1 | 0x4, rate_value)
        v.require(k.SetInformationJobObject(job, 15, ct.byref(rate), ct.sizeof(rate)))
        info = STARTUPINFOEX()
        info.StartupInfo.cb = ct.sizeof(info)
        info.StartupInfo.dwFlags = 1 | 0x100
        info.StartupInfo.wShowWindow = 0
        handles = (wt.HANDLE * 2)(msvcrt.get_osfhandle(stdin.fileno()), msvcrt.get_osfhandle(console.fileno()))
        info.StartupInfo.hStdInput = handles[0]
        info.StartupInfo.hStdOutput = info.StartupInfo.hStdError = handles[1]
        size = ct.c_size_t()
        k.InitializeProcThreadAttributeList(None, 1, 0, ct.byref(size))
        if not size.value:
            raise ct.WinError(ct.get_last_error())
        attributes = ct.create_string_buffer(size.value)
        v.require(k.InitializeProcThreadAttributeList(attributes, 1, 0, ct.byref(size)))
        initialized = True
        info.lpAttributeList = ct.cast(attributes, ct.c_void_p)
        environment = dict(os.environ, GAUSS_EXEDIR=str(Path(argv[0]).parent),
                           GAUSS_SCRDIR=str(cwd / 'scratch'), OMP_NUM_THREADS=str(cores),
                           NCPUS=str(cores), OMP_THREAD_LIMIT=str(cores),
                           MKL_NUM_THREADS=str(cores), OPENBLAS_NUM_THREADS='1')
        env = ct.create_unicode_buffer('\0'.join(f'{a}={b}' for a, b in sorted(environment.items(),
                                      key=lambda item: item[0].upper())) + '\0\0')
        command = ct.create_unicode_buffer(subprocess.list2cmdline([str(x) for x in argv]))
        # Restricted handle list avoids one worker inheriting another's files.
        with _LAUNCH_LOCK:
            for handle in handles:
                os.set_handle_inheritable(handle, True)
            try:
                v.require(k.UpdateProcThreadAttribute(attributes, 0, 0x20002, handles, ct.sizeof(handles), None, None))
                v.require(k.CreateProcessW(str(argv[0]), command, None, None, True,
                          0x4 | 0x400 | 0x8 | 0x80000, env, str(cwd), ct.byref(info), ct.byref(proc)))
            finally:
                for handle in handles:
                    os.set_handle_inheritable(handle, False)
        try:
            v.require(k.AssignProcessToJobObject(job, proc.hProcess))
        except BaseException:
            v.require(k.TerminateProcess(proc.hProcess, 125))
            if k.WaitForSingleObject(proc.hProcess, 5000) != 0:
                raise RuntimeError('Unassigned owned process did not confirm exit')
            raise
        on_started({'pid': proc.dwProcessId, 'started_epoch': started, 'cores': cores,
                    'cpu_rate_hard_cap': rate_value, 'memory_limit_gib': memory_gib,
                    'console_policy': 'detached; system console infrastructure remains responsive'})
        if k.ResumeThread(proc.hThread) == 0xFFFFFFFF:
            raise ct.WinError(ct.get_last_error())
        accounting, last_tick = v.ACCOUNTING(), 0.0
        def emit(state=None, diagnostic=None):
            current = time.monotonic()
            v.require(k.QueryInformationJobObject(job, 1, ct.byref(accounting), ct.sizeof(accounting), None))
            on_tick({'elapsed_seconds': current - begin,
                'active_elapsed_seconds': clock.active(current), 'ram_paused_seconds': clock.paused(current),
                'active_processes': accounting.ActiveProcesses,
                'execution_state': state or ('ram_paused' if pause.held else pause.state if pause.state == 'pause_failed' else 'running'),
                'pause_is_disk_checkpoint': False, 'suspended_thread_handles': len(pause.handles),
                'pause_error': pause_error, 'pause_diagnostic': diagnostic or pause.last_diagnostic,
                'active_timeout_count': timeout_count,
                'tree_cpu_seconds': (accounting.TotalUserTime + accounting.TotalKernelTime) / 1e7})
        while True:
            v.require(k.QueryInformationJobObject(job, 1, ct.byref(accounting), ct.sizeof(accounting), None))
            if not accounting.ActiveProcesses:
                break
            now = time.monotonic()
            # Exhausted active time parks the computation when the caller can
            # persist a hold command. Utilities retain their bounded timeout.
            if reason is None and clock.active(now) >= next_timeout and not timeout_notified:
                timeout_notified = True
                timeout_count += 1
                if on_timeout is not None:
                    on_timeout()
                else:
                    reason = 'timeout'
            if reason is None and cancel():
                reason = 'operator_interrupt'
            if reason is not None:
                if termination_started is None:
                    v.require(k.TerminateJobObject(job, 123 if reason == 'operator_interrupt' else 124))
                    termination_started = now
                    emit('terminating')
                elif now - termination_started >= 5:
                    raise RuntimeError('Owned Job termination did not complete before its deadline')
            else:
                wanted = bool(hold())
                if not wanted:
                    hold_failed = False
                    if previous_wanted and timeout_notified:
                        # An explicit release grants another bounded active interval.
                        next_timeout = clock.active(now) + timeout
                        timeout_notified = False
                previous_hold = pause.held
                if wanted and not hold_failed and (not pause.held or now - last_pause_scan > 1):
                    try:
                        pause.hold(refresh=pause.held, cancel=lambda: cancel() or not hold(),
                                   on_progress=lambda detail: emit('pausing', detail))
                        last_pause_scan = time.monotonic()
                        pause_error = None
                    except Exception as error:
                        # hold() rolls back partial suspensions. Never silently
                        # turn a failed RAM pause into a kill of the calculation.
                        if pause.held:
                            raise
                        hold_failed = True
                        pause_error = str(error)
                        pause_failures.append(dict(pause.last_diagnostic, epoch=time.time()))
                        with (cwd / 'pause-failures.jsonl').open('a', encoding='utf-8') as log:
                            log.write(json.dumps(pause_failures[-1]) + '\n')
                        on_hold_error(pause_error)
                elif not wanted and pause.held:
                    pause.resume()
                elif not wanted:
                    pause.state = 'running'
                clock.set_held(pause.held, time.monotonic())
                previous_wanted = wanted
                if previous_hold != pause.held:
                    last_tick = 0.0
            now = time.monotonic()
            if now - last_tick >= .5:
                emit('terminating' if reason else None)
                last_tick = now
            time.sleep(0.1)
        exit_code = wt.DWORD()
        v.require(k.GetExitCodeProcess(proc.hProcess, ct.byref(exit_code)))
        v.require(k.QueryInformationJobObject(job, 9, ct.byref(limits), ct.sizeof(limits), None))
        return {'pid': proc.dwProcessId, 'argv': [str(x) for x in argv], 'started_epoch': started,
                'finished_epoch': time.time(), 'wall_seconds': time.monotonic() - begin,
                'active_wall_seconds': clock.active(time.monotonic()),
                'ram_paused_seconds': clock.paused(time.monotonic()),
                'tree_cpu_seconds': (accounting.TotalUserTime + accounting.TotalKernelTime) / 1e7,
                'peak_tree_commit_bytes': limits.PeakJobMemoryUsed, 'cores': cores,
                'memory_limit_gib': memory_gib, 'cpu_rate_hard_cap': rate_value,
                 'tree_process_count': accounting.TotalProcesses, 'active_processes_at_return': accounting.ActiveProcesses,
                 'active_timeout_count': timeout_count, 'pause_failures': pause_failures,
                  'console_policy': 'detached; System32 conhost excluded from compute suspension',
                'exit_code': exit_code.value, 'interrupted': reason is not None, 'stop_reason': reason}
    finally:
        try:
            _drain_owned_tree(k, v, job)
        finally:
            # KILL_ON_JOB_CLOSE remains the final ownership backstop on any error.
            k.CloseHandle(job)
            if pause is not None:
                for handle in pause.handles.values():
                    k.CloseHandle(handle)
                pause.handles.clear()
            if proc.hThread:
                k.CloseHandle(proc.hThread)
            if proc.hProcess:
                k.CloseHandle(proc.hProcess)
            if initialized:
                k.DeleteProcThreadAttributeList(attributes)
            if stdin is not None:
                stdin.close()
            if console is not None:
                console.close()
