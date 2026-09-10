"""The scientific root gate must stop the medium launch before any configuration/run mutation."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import run_steady_acceleration as runner
import analyze_steady_acceleration as analysis


class SteadyGuardTests(unittest.TestCase):
    def test_failed_coarse_equivalence_blocks_medium_before_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run = root/'runs/steady-acceleration-c2c'
            run.mkdir(parents=True)
            ref = root/'runs/scientific-closure/128x32/attempt-001'
            ref.mkdir(parents=True)
            source = ref/'config.json'
            source.write_text('{}')
            (run/'preset.json').write_text(json.dumps({
                'coarse_source_sha256':{'config.json':hashlib.sha256(source.read_bytes()).hexdigest()}}))
            with patch.object(runner,'ROOT',root), patch.object(runner,'RUN',run), \
                    patch.object(sys,'argv',['run_steady_acceleration.py','medium']), \
                    patch.object(analysis,'coarse_evidence',return_value={'root_equivalence_pass':False}), \
                    patch.object(runner.subprocess,'run') as launch:
                with self.assertRaisesRegex(RuntimeError,'medium is prohibited'):
                    runner.main()
                launch.assert_not_called()
            self.assertFalse((run/'medium-config.json').exists())
            self.assertFalse((run/'medium-256x64').exists())

    def test_reference_change_stops_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            run=root/'runs/steady-acceleration-c2c';run.mkdir(parents=True)
            ref=root/'runs/scientific-closure/128x32/attempt-001';ref.mkdir(parents=True)
            (ref/'config.json').write_text('{}')
            (run/'preset.json').write_text(json.dumps({'coarse_source_sha256':{'config.json':'bad'}}))
            with patch.object(runner,'ROOT',root),patch.object(runner,'RUN',run), \
                    patch.object(sys,'argv',['run_steady_acceleration.py','coarse']), \
                    patch.object(runner.subprocess,'run') as launch:
                with self.assertRaisesRegex(ValueError,'reference changed'):
                    runner.main()
                launch.assert_not_called()


if __name__=='__main__':
    unittest.main()
