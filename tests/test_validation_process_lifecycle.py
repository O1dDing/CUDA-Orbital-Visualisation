"""Real Windows process checks for lifecycle timing; evidence is retained."""
from __future__ import annotations
import json
import os
from pathlib import Path
import sys
import tempfile
import time
import unittest

from validation_process import atomic_json, run_tree


@unittest.skipUnless(os.name == 'nt', 'Windows Job Objects are required')
class ValidationProcessLifecycleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        specified = os.environ.get('COV_SUPERVISOR_TEST_ROOT')
        cls.root = Path(specified).resolve() if specified else Path(tempfile.mkdtemp(prefix='cov-process-lifecycle-')).resolve()
        cls.root.mkdir(parents=True, exist_ok=True)
        cls.program = cls.root / 'probe child.py'
        cls.program.write_text('''import json, os, pathlib, sys, time
value = {'pid': os.getpid(), 'input': sys.stdin.read(), 'argument': sys.argv[2],
         'marker': os.environ['COV_LIFECYCLE_TEST_MARKER']}
pathlib.Path(sys.argv[1]).write_text(json.dumps(value, ensure_ascii=False), encoding='utf-8')
print('LIFECYCLE_STDOUT_OK', flush=True)
if sys.argv[3] == 'sleep':
    time.sleep(60)
elif sys.argv[3] == 'fail':
    raise SystemExit(7)
''', encoding='utf-8')

    def run_probe(self, name, mode='normal', on_phase=None, on_started=None, timeout=120):
        directory = self.root / name
        directory.mkdir(exist_ok=False)
        source = directory / 'input.txt'
        source.write_text('输入与标准句柄\nsecond line\n', encoding='utf-8')
        output = directory / '输出 with space.json'
        argument = 'quoted "value" $() ; & `'
        env = dict(os.environ, PYTHONIOENCODING='utf-8', OMP_NUM_THREADS='1',
                   OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', COV_LIFECYCLE_TEST_MARKER=name)
        result = run_tree([sys.executable, '-X', 'utf8', str(self.program), str(output), argument, mode],
                          directory, env, 1, 4, timeout, on_started, affinity=False,
                          stdin_path=source, on_phase=on_phase)
        atomic_json(directory / 'process.json', result)
        return result, output, argument

    def assert_complete_lifecycle(self, result):
        events = result['lifecycle_events']
        names = [event['phase'] for event in events]
        self.assertEqual(names[0], 'supervisor_enter')
        self.assertEqual(names[-1], 'cleanup_complete')
        sequence = ['process_creation_enter', 'process_created_suspended', 'job_assignment_enter',
                    'job_assigned', 'resume_thread_enter', 'thread_resumed', 'tree_monitor_enter',
                    'tree_empty', 'final_accounting_collected', 'cleanup_enter', 'cleanup_complete']
        self.assertEqual([names.index(name) for name in sequence], sorted(names.index(name) for name in sequence))
        for index, event in enumerate(events):
            self.assertEqual(event['sequence'], index)
            self.assertGreaterEqual(event['observer_wall_seconds'], 0)
            self.assertGreaterEqual(event['after_observer_elapsed_seconds'], event['elapsed_seconds'])
            if index:
                self.assertGreaterEqual(event['elapsed_seconds'], events[index-1]['after_observer_elapsed_seconds'])
        self.assertEqual(result['core_count'], 1)
        self.assertEqual(result['memory_limit_gib'], 4)
        self.assertIsNone(result['cpu_mask'])
        self.assertGreater(result['cpu_rate_hard_cap'], 0)

    def test_real_handles_argument_roundtrip_and_observer_overhead(self):
        observed = []
        def observe(event):
            observed.append(event)
            if event['phase'] == 'process_creation_enter':
                time.sleep(.03)
            if event['phase'] == 'job_assigned':
                raise RuntimeError('intentional telemetry failure')
        def started(info):
            self.assertGreater(info['pid'], 0)
            time.sleep(.04)
        result, output, argument = self.run_probe('normal 句柄', on_phase=observe, on_started=started)
        self.assertEqual(result['exit_code'], 0)
        self.assertFalse(result['timed_out'])
        self.assert_complete_lifecycle(result)
        child = json.loads(output.read_text(encoding='utf-8'))
        self.assertEqual(child['input'], '输入与标准句柄\nsecond line\n')
        self.assertEqual(child['argument'], argument)
        self.assertEqual(child['marker'], 'normal 句柄')
        self.assertIn('LIFECYCLE_STDOUT_OK', output.parent.joinpath('launcher.log').read_text(encoding='utf-8'))
        phases = {event['phase']: event for event in result['lifecycle_events']}
        self.assertGreaterEqual(phases['process_creation_enter']['observer_wall_seconds'], .02)
        self.assertGreaterEqual(phases['process_created_suspended']['elapsed_seconds'],
                                phases['process_creation_enter']['after_observer_elapsed_seconds'])
        self.assertGreaterEqual(phases['started_callback_returned']['elapsed_seconds']-
                                phases['started_callback_enter']['after_observer_elapsed_seconds'], .03)
        self.assertEqual(phases['job_assigned']['observer_error'], 'RuntimeError: intentional telemetry failure')
        self.assertEqual(len(observed), len(result['lifecycle_events']))

    def test_nonzero_child_exit_retains_failure_and_cleanup_with_default_observer(self):
        result, output, _ = self.run_probe('exit seven', mode='fail')
        self.assertEqual(result['exit_code'], 7)
        self.assertFalse(result['timed_out'])
        self.assertTrue(output.exists())
        self.assert_complete_lifecycle(result)
        self.assertTrue(all('observer_error' not in event for event in result['lifecycle_events']))
        self.assertIn('started_callback_absent', [event['phase'] for event in result['lifecycle_events']])

    def test_actual_execution_timeout_records_deadline_and_empty_job_before_cleanup(self):
        result, output, _ = self.run_probe('execution deadline', mode='sleep', timeout=2)
        self.assertTrue(result['timed_out'])
        self.assertEqual(result['exit_code'], 124)
        self.assertTrue(output.exists(), 'The child must have actually executed before this timeout check is credited')
        self.assert_complete_lifecycle(result)
        names = [event['phase'] for event in result['lifecycle_events']]
        self.assertEqual(names.count('execution_deadline_observed'), 1)
        self.assertLess(names.index('execution_deadline_observed'), names.index('tree_empty'))


if __name__ == '__main__':
    unittest.main()
