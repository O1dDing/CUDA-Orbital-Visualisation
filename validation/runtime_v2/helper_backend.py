"""Desktop Helper adapter; construction and status refresh never dispatch work."""
from pathlib import Path
import json
import os
import re
import subprocess
import sys
import time

from state_store import read_json

HELPER_VERSION = '2.2.1'
CONTROL_COMMANDS = frozenset(('hold', 'pause', 'resume', 'interrupt', 'shutdown'))


class HelperBackend:
    def __init__(self, home, *, runner=subprocess.run, launcher=subprocess.Popen):
        self.home = Path(home).resolve()
        self.pointer = read_json(self.home / 'runtime/active.json')
        self.bundle = Path(self.pointer['bundle']).resolve()
        if not self.bundle.is_relative_to(self.home / 'runtime/releases'):
            raise ValueError('Active runtime must belong to this deployment')
        self.config_path = self.home / 'runtime/config.json'
        if Path(self.pointer['config']).resolve() != self.config_path:
            raise ValueError('The Helper requires the shared deployment configuration')
        self.config = read_json(self.config_path)
        directory = self.config['runtime_directory']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,95}', directory):
            raise ValueError('Invalid runtime directory')
        self.work = Path(self.config['work_root']).resolve()
        self.runtime_root = self.work / directory
        self.runtime = self.bundle / 'validation/runtime_v2/resume.py'
        if not self.runtime.is_file():
            raise FileNotFoundError(self.runtime)
        executable = Path(sys.executable)
        self.python = executable.with_name('python.exe') if executable.name.lower() == 'pythonw.exe' else executable
        self.runner, self.launcher = runner, launcher
        # This permission exists only in this window's memory. It is never saved.
        self.armed = False

    def argv(self, command, *extra):
        return [str(self.python), '-B', '-X', 'utf8', str(self.runtime),
                '--config', str(self.config_path), command, *extra]

    def _short_command(self, command):
        result = self.runner(self.argv(command), capture_output=True, encoding='utf-8',
                             timeout=30, cwd=str(self.home),
                             creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        if result.returncode:
            raise RuntimeError((result.stderr or result.stdout).strip())
        return result.stdout

    def status(self):
        raw = self._short_command('status')
        status, end = json.JSONDecoder().raw_decode(raw.lstrip())
        status['progress_text'] = raw.lstrip()[end:].strip()
        status['binding'] = read_json(self.runtime_root / 'binding.json', {})
        return status

    def check(self):
        value = json.loads(self._short_command('check'))
        if value.get('calculation_started') is not False:
            raise RuntimeError('Deployment check did not confirm a non-computing result')
        return value

    def arm(self, enabled):
        self.armed = bool(enabled)

    def _require_armed(self):
        if not self.armed:
            raise RuntimeError('Enable manual calculation controls in this window first')

    def control(self, command):
        self._require_armed()
        if command not in CONTROL_COMMANDS:
            raise ValueError('Unsupported control command')
        state = self.status()
        if not state.get('work_lease_held'):
            raise RuntimeError('No coordinator is running; no control command was written')
        return self._short_command(command)

    def start(self, limit, workers, cases=''):
        self._require_armed()
        limit, workers = int(limit), int(workers)
        if not 1 <= limit <= 273 or workers not in (1, 2):
            raise ValueError('Case limit must be 1–273; worker slots must be 1 or 2')
        selected = list(dict.fromkeys(c.upper() for c in re.split(r'[\s,;]+', cases.strip()) if c))
        if any(not re.fullmatch(r'OLD-\d{3}', c) or not 1 <= int(c[4:]) <= 273 for c in selected):
            raise ValueError('Use case IDs from OLD-001 through OLD-273')
        if self.status().get('work_lease_held'):
            raise RuntimeError('A coordinator already owns this calculation directory')
        extra = ['--limit', str(limit), '--workers', str(workers)]
        for case in selected:
            extra.extend(['--case', case])
        logs = self.home / 'runtime/helper-logs'
        logs.mkdir(exist_ok=True)
        log = logs / ('coordinator-' + str(time.time_ns()) + '.log')
        with log.open('xb') as stream:
            process = self.launcher(self.argv('run', *extra), cwd=str(self.home),
                                    stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT,
                                    creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0) |
                                                  getattr(subprocess, 'CREATE_NEW_PROCESS_GROUP', 0))
        return {'pid': process.pid, 'log': str(log)}
