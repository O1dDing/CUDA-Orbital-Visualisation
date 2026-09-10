"""Reversible RAM-only pause of threads in ONE owned Windows Job Object.

SuspendThread is debugger-style control, not a checkpoint transaction. Never
copy/convert Gaussian files or acquire a target's locks while it is held. Kernel
I/O may still finish. Keep every successful suspend paired with one resume;
never reset somebody else's suspend count. No process-name matching or Nt APIs.
"""
from __future__ import annotations

import ctypes as ct
from ctypes import wintypes as wt
import time


class THREADENTRY32(ct.Structure):
    _fields_ = [(name, wt.DWORD) for name in ('dwSize', 'cntUsage', 'th32ThreadID', 'th32OwnerProcessID')] + [
        ('tpBasePri', wt.LONG), ('tpDeltaPri', wt.LONG), ('dwFlags', wt.DWORD)]


class JobPause:
    def __init__(self, kernel, job):
        self.k, self.job = kernel, job
        self.handles = {}  # (tid, creation FILETIME) -> HANDLE; exactly one owned increment
        self.held = False
        signatures = {
            'OpenProcess': ([wt.DWORD, wt.BOOL, wt.DWORD], wt.HANDLE),
            'OpenThread': ([wt.DWORD, wt.BOOL, wt.DWORD], wt.HANDLE),
            'IsProcessInJob': ([wt.HANDLE, wt.HANDLE, ct.POINTER(wt.BOOL)], wt.BOOL),
            'GetProcessIdOfThread': ([wt.HANDLE], wt.DWORD),
            'GetThreadTimes': ([wt.HANDLE] + [ct.POINTER(wt.FILETIME)] * 4, wt.BOOL),
            'SuspendThread': ([wt.HANDLE], wt.DWORD),
            'ResumeThread': ([wt.HANDLE], wt.DWORD),
            'WaitForSingleObject': ([wt.HANDLE, wt.DWORD], wt.DWORD),
            'CreateToolhelp32Snapshot': ([wt.DWORD, wt.DWORD], wt.HANDLE),
            'Thread32First': ([wt.HANDLE, ct.POINTER(THREADENTRY32)], wt.BOOL),
            'Thread32Next': ([wt.HANDLE, ct.POINTER(THREADENTRY32)], wt.BOOL),
            'CloseHandle': ([wt.HANDLE], wt.BOOL),
        }
        for name, (args, result) in signatures.items():
            fun = getattr(kernel, name)
            fun.argtypes, fun.restype = args, result

    def pids(self):
        capacity = 16
        while capacity <= 65536:
            class IDS(ct.Structure):
                _fields_ = [('assigned', wt.DWORD), ('listed', wt.DWORD), ('pids', ct.c_size_t * capacity)]
            ids = IDS()
            ok = self.k.QueryInformationJobObject(self.job, 3, ct.byref(ids), ct.sizeof(ids), None)
            if ok and ids.listed >= ids.assigned:
                return set(ids.pids[:ids.listed])
            if not ok and ct.get_last_error() != 234:  # ERROR_MORE_DATA
                raise ct.WinError(ct.get_last_error())
            capacity = max(capacity * 2, int(ids.assigned) + 8)
        raise RuntimeError('Job process list too large to pause safely')

    def threads(self, pids):
        snapshot = self.k.CreateToolhelp32Snapshot(0x4, 0)  # TH32CS_SNAPTHREAD
        if snapshot == ct.c_void_p(-1).value:
            raise ct.WinError(ct.get_last_error())
        found = []
        try:
            entry = THREADENTRY32()
            entry.dwSize = ct.sizeof(entry)
            ok = self.k.Thread32First(snapshot, ct.byref(entry))
            while ok:
                if entry.th32OwnerProcessID in pids:
                    found.append((int(entry.th32OwnerProcessID), int(entry.th32ThreadID)))
                entry.dwSize = ct.sizeof(entry)
                ok = self.k.Thread32Next(snapshot, ct.byref(entry))
            if ct.get_last_error() != 18:  # ERROR_NO_MORE_FILES
                raise ct.WinError(ct.get_last_error())
            return found
        finally:
            self.k.CloseHandle(snapshot)

    def hold(self, timeout=5.0, refresh=False):
        if self.held and not refresh:
            return len(self.handles)
        begin, stable, previous = time.monotonic(), 0, None
        try:
            while time.monotonic() - begin < timeout:
                pids = self.pids()
                live = set()
                for pid, tid in self.threads(pids):
                    # Pin the thread object before checking its actual owner/Job.
                    handle = self.k.OpenThread(0x2 | 0x40 | 0x100000, False, tid)
                    if not handle:
                        if ct.get_last_error() == 87:  # vanished between snapshot and open
                            continue
                        raise ct.WinError(ct.get_last_error())
                    retained = False
                    try:
                        if self.k.GetProcessIdOfThread(handle) != pid:
                            continue
                        process = self.k.OpenProcess(0x1000, False, pid)
                        if not process:
                            if self.k.WaitForSingleObject(handle, 0) == 0:
                                continue
                            raise ct.WinError(ct.get_last_error())
                        try:
                            member = wt.BOOL()
                            if not self.k.IsProcessInJob(process, self.job, ct.byref(member)):
                                raise ct.WinError(ct.get_last_error())
                            if not member.value:
                                continue
                        finally:
                            self.k.CloseHandle(process)
                        times = [wt.FILETIME() for _ in range(4)]
                        if not self.k.GetThreadTimes(handle, *(ct.byref(t) for t in times)):
                            raise ct.WinError(ct.get_last_error())
                        created = (times[0].dwHighDateTime << 32) | times[0].dwLowDateTime
                        key = (tid, created)
                        if self.k.WaitForSingleObject(handle, 0) == 0:
                            continue
                        live.add(key)
                        if key not in self.handles:
                            if self.k.SuspendThread(handle) == 0xFFFFFFFF:
                                if self.k.WaitForSingleObject(handle, 0) == 0:
                                    live.discard(key)
                                    continue
                                raise ct.WinError(ct.get_last_error())
                            self.handles[key] = handle
                            retained = True
                    finally:
                        if not retained:
                            self.k.CloseHandle(handle)
                observed = (frozenset(pids), frozenset(live))
                stable = stable + 1 if observed == previous else 0
                previous = observed
                if stable >= 2 and (not pids or live):
                    self.held = True
                    return len(self.handles)
                time.sleep(.04)
            raise RuntimeError('Could not stabilize owned thread tree within pause deadline')
        except BaseException:
            self.resume()  # rollback every partial suspension, or raise if rollback fails
            raise

    def resume(self):
        errors = []
        for key, handle in list(self.handles.items()):
            exited = self.k.WaitForSingleObject(handle, 0) == 0
            if not exited and self.k.ResumeThread(handle) == 0xFFFFFFFF:
                if self.k.WaitForSingleObject(handle, 0) != 0:
                    errors.append(key)
                    continue
            self.k.CloseHandle(handle)
            del self.handles[key]
        self.held = bool(self.handles)
        if errors:
            raise RuntimeError('Could not undo owned suspend counts: ' + repr(errors))

    def close_after_exit(self):
        # Used only after the owning supervisor confirmed zero active processes.
        for handle in self.handles.values():
            self.k.CloseHandle(handle)
        self.handles.clear()
        self.held = False
