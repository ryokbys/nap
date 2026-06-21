import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

__version__ = '260621'


class TestPMD(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parents[2]
        cls.py_run_dir = root / 'test' / 'tmp_pmd' / 'py_run'

    def _prepare_workdir(self):
        tmpdir = tempfile.TemporaryDirectory()
        workdir = Path(tmpdir.name)
        for name in ('in.pmd', 'in.params.Coulomb', 'in.params.Morse', 'pmdini'):
            shutil.copy2(self.py_run_dir / name, workdir / name)
        return tmpdir, workdir

    def _run_documented_flow(self, *, pass_system_in_ctor):
        tmpdir, workdir = self._prepare_workdir()
        try:
            code = f"""
import json
import numpy as np
import nappy
from nappy.pmd import PMD

nsys = nappy.io.read('pmdini')
before = nsys.get_scaled_positions().copy()
if {pass_system_in_ctor!r}:
    pmd = PMD(nsys)
    pmd.load_inpmd()
else:
    pmd = PMD()
    pmd.load_inpmd()
    pmd.set_system(nsys)
pmd.run()
out = pmd.get_system()
print(json.dumps({{
    'disp': float(np.linalg.norm(out.get_scaled_positions() - before)),
    'epot': out.get_potential_energy(),
}}))
"""
            proc = subprocess.run(
                [sys.executable, '-c', code],
                cwd=workdir,
                text=True,
                capture_output=True,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)
            result = json.loads(proc.stdout.splitlines()[-1])
            return result
        finally:
            tmpdir.cleanup()

    def test_run_uses_loaded_inpmd_and_system_state(self):
        result = self._run_documented_flow(pass_system_in_ctor=True)
        self.assertGreater(result['disp'], 1.0e-6)

    def test_set_system_after_load_inpmd_keeps_specorder(self):
        ctor_result = self._run_documented_flow(pass_system_in_ctor=True)
        set_result = self._run_documented_flow(pass_system_in_ctor=False)
        self.assertGreater(set_result['disp'], 1.0e-6)
        self.assertAlmostEqual(set_result['disp'], ctor_result['disp'], places=10)
        self.assertAlmostEqual(set_result['epot'], ctor_result['epot'], places=10)


if __name__ == '__main__':
    unittest.main()
