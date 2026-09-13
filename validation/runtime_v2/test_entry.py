import tempfile
from pathlib import Path
import unittest
import os
from unittest.mock import patch
from entry import translate
from install_runtime import stub, activate, rollback
import resume
from test_runtime import FakeData, FakeBackend


class EntryTests(unittest.TestCase):
    def test_legacy_flags_preserve_selection_limit_and_worker_count(self):
        self.assertEqual(translate(['--legacy']), ['menu'])
        self.assertEqual(translate(['--legacy', '--execute', '--limit', '2', '--workers', '1', '--case', 'OLD-093']),
            ['run', '--limit', '2', '--workers', '1', '--case', 'OLD-093'])
        self.assertEqual(translate(['status']), ['status'])

    def test_modern_and_legacy_stubs_parse_with_space_and_unicode_paths(self):
        for legacy in (True, False):
            content = stub(Path('C:/COV 空格/Resume'), 'resume.py', legacy)
            compile(content, '<compatible entry>', 'exec')

    def test_activation_and_rollback_preserve_original_entries_and_results(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            home, fast = root / 'Resume', root / 'FastPause-2.1'
            data = FakeData(root)
            config = {'work_root': str(root / 'work')}
            engine = resume.Engine(config, data, FakeBackend()); engine.bind()
            resume.atomic(engine.root / 'import-summary.json', {'test_fixture': True})
            proof = root / 'proof.txt'; proof.write_text('fixture native evidence')
            resume.atomic(engine.root / 'native-capabilities.json', {'runtime_identity': engine.identity,
                'binaries': data.binaries, 'capabilities': {'ram_pause': {'passed': True,
                    'methods': ['RPBE1PBE', 'UPBE1PBE']}},
                'evidence': [{'path': str(proof), 'sha256': resume.sha(proof)}]})
            config_file = root / 'config.json'; resume.atomic(config_file, config)
            paths = [home / 'cov_resume_menu.py', home / 'README-先读我.md'] + [
                fast / 'validation/runtime_v2' / f for f in ('resume.py', 'config.json', 'README.md', 'FAST_PAUSE.md', 'WORK_PROMPT.md')]
            for path in paths:
                path.parent.mkdir(parents=True, exist_ok=True); path.write_text('original ' + path.name)
            before = {p: p.read_bytes() for p in paths}
            prepared = {'bundle_identity': 'fixture', 'bundle': str(root / 'bundle'), 'config': str(config_file)}
            with patch.dict(os.environ, {'LOCALAPPDATA': str(root / 'local')}), \
                    patch.object(resume, 'FrozenData', return_value=data), \
                    patch('windows_job.existing_gaussian', return_value=[]):
                result = activate(home, fast, prepared)
                self.assertIn('active.json', (home / 'cov_resume_menu.py').read_text())
                receipt = Path(result['rollback_receipt'])
                (home / 'cov_resume_menu.py').write_text('later edit')
                with self.assertRaisesRegex(ValueError, 'Entry changed'):
                    rollback(home, receipt)
                (home / 'cov_resume_menu.py').write_bytes(stub(home, 'resume.py', True).encode('utf-8'))
                restored = rollback(home, receipt)
                self.assertEqual(restored['restored_files'], len(paths))
                self.assertEqual(before, {p: p.read_bytes() for p in paths})
                self.assertEqual(proof.read_text(), 'fixture native evidence')


if __name__ == '__main__':
    unittest.main(verbosity=2)
