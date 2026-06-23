"""
Unit tests for NAPSystem.make_pair_list() – scipy fallback path.

The Fortran pairlist module (nappy.pmd.pairlist) is disabled by patching
sys.modules so only the scipy.spatial.cKDTree code path is exercised.
The cross-validation test re-enables the Fortran module and compares both
results to confirm parity.
"""
import sys
import unittest
import numpy as np

__version__ = '260623'


def _make_fcc(a=4.0, n=3):
    """FCC n×n×n supercell; actual cell length = a*n Å."""
    # Import here so the module-level Fortran patch doesn't interfere.
    from nappy.napsys import NAPSystem
    nsys = NAPSystem()
    nsys.set_lattice(
        a * n,
        np.array([1., 0., 0.]),
        np.array([0., 1., 0.]),
        np.array([0., 0., 1.]))
    basis = [[0., 0., 0.], [.5, .5, 0.], [.5, 0., .5], [0., .5, .5]]
    poss = []
    for i1 in range(n):
        for i2 in range(n):
            for i3 in range(n):
                for b in basis:
                    poss.append([(b[0]+i1)/n, (b[1]+i2)/n, (b[2]+i3)/n])
    nsys.add_atoms(['Al'] * len(poss), np.array(poss))
    return nsys


class TestMakePairListScipy(unittest.TestCase):
    """scipy cKDTree fallback for NAPSystem.make_pair_list()."""

    def setUp(self):
        """Disable the Fortran pairlist module to force the scipy path."""
        self._orig = sys.modules.get('nappy.pmd.pairlist')
        sys.modules['nappy.pmd.pairlist'] = None

    def tearDown(self):
        """Restore the Fortran pairlist module."""
        if self._orig is None:
            sys.modules.pop('nappy.pmd.pairlist', None)
        else:
            sys.modules['nappy.pmd.pairlist'] = self._orig

    # ── helpers ──────────────────────────────────────────────────────────────

    def _check_symmetric(self, nsys):
        natm = nsys.num_atoms()
        for ia in range(natm):
            for ja in nsys.atoms['neighbors'][ia]:
                self.assertIn(ia, nsys.atoms['neighbors'][ja],
                              f'asymmetric: {ia} lists {ja} but not vice versa')

    def _check_no_duplicates(self, nsys):
        natm = nsys.num_atoms()
        for ia in range(natm):
            nb = nsys.atoms['neighbors'][ia]
            self.assertEqual(len(nb), len(set(nb)),
                             f'duplicate neighbor at atom {ia}')

    # ── tests ─────────────────────────────────────────────────────────────────

    def test_fcc_nn_count(self):
        """FCC 3×3×3: every atom has exactly 12 NN within 3.2 Å."""
        nsys = _make_fcc(a=4.0, n=3)
        nsys.make_pair_list(rcut=3.2)
        counts = [len(nsys.atoms['neighbors'][i]) for i in range(nsys.num_atoms())]
        self.assertTrue(all(c == 12 for c in counts),
                        f'NN counts: min={min(counts)} max={max(counts)}')

    def test_symmetric(self):
        """Neighbor list must be symmetric: ia in nbrs[ja] iff ja in nbrs[ia]."""
        nsys = _make_fcc(a=4.0, n=3)
        nsys.make_pair_list(rcut=3.2)
        self._check_symmetric(nsys)

    def test_no_duplicates(self):
        """No atom index should appear twice in another atom's neighbor list."""
        nsys = _make_fcc(a=4.0, n=3)
        nsys.make_pair_list(rcut=3.2)
        self._check_no_duplicates(nsys)

    def test_multiple_images_deduplicated(self):
        """4-atom FCC cell (cell < 2*rcut): 4 images per NN pair, stored once."""
        from nappy.napsys import NAPSystem
        a = 4.0
        nsys = NAPSystem()
        nsys.set_lattice(a,
                         np.array([1., 0., 0.]),
                         np.array([0., 1., 0.]),
                         np.array([0., 0., 1.]))
        nsys.add_atoms(['Al']*4,
                       np.array([[0.,0.,0.],[.5,.5,0.],[.5,0.,.5],[0.,.5,.5]]))
        nsys.make_pair_list(rcut=2.9)
        counts = [len(nsys.atoms['neighbors'][i]) for i in range(4)]
        # Each atom has 3 unique neighbor atoms (12 NN positions / 4 images each)
        self.assertEqual(counts, [3, 3, 3, 3],
                         f'Expected [3,3,3,3], got {counts}')
        self._check_symmetric(nsys)
        self._check_no_duplicates(nsys)

    def test_non_orthogonal_cell(self):
        """Sheared FCC 3×3×3: NN count must remain 12 per atom."""
        nsys = _make_fcc(a=4.0, n=3)
        # Apply a 0.5 Å shear to a2
        nsys.a2[0] += 0.5 / (4.0 * 3)
        nsys.make_pair_list(rcut=3.2)
        counts = [len(nsys.atoms['neighbors'][i]) for i in range(nsys.num_atoms())]
        self.assertTrue(all(c == 12 for c in counts),
                        f'Sheared FCC NN counts: {set(counts)}')
        self._check_symmetric(nsys)

    def test_pairwise_rcuts(self):
        """Species-specific cutoffs: Si-O pairs excluded when cutoff is tiny."""
        from nappy.napsys import NAPSystem
        n = 3
        a = 4.0
        nsys = NAPSystem()
        nsys.set_lattice(a * n,
                         np.array([1., 0., 0.]),
                         np.array([0., 1., 0.]),
                         np.array([0., 0., 1.]))
        basis = [[0.,0.,0.],[.5,.5,0.],[.5,0.,.5],[0.,.5,.5]]
        poss, sps = [], []
        for i1 in range(n):
            for i2 in range(n):
                for i3 in range(n):
                    for k, b in enumerate(basis):
                        poss.append([(b[0]+i1)/n,(b[1]+i2)/n,(b[2]+i3)/n])
                        sps.append('Si' if k % 2 == 0 else 'O')
        nsys.add_atoms(sps, np.array(poss))
        nsys.make_pair_list(rcut=3.2,
                            rcuts={('Si','Si'):3.2, ('O','O'):3.2, ('Si','O'):0.1})
        for ia in range(nsys.num_atoms()):
            sp_ia = nsys.specorder[nsys.atoms.sid[ia]-1]
            for ja in nsys.atoms['neighbors'][ia]:
                sp_ja = nsys.specorder[nsys.atoms.sid[ja]-1]
                self.assertEqual(sp_ia, sp_ja,
                                 f'Cross-species pair {ia}({sp_ia})-{ja}({sp_ja})')


class TestMakePairListCrossValidation(unittest.TestCase):
    """Cross-validate scipy path against the Fortran path."""

    def test_matches_fortran(self):
        """scipy and Fortran paths must produce identical sorted neighbor lists."""
        # Check if the Fortran module is actually compiled and usable
        try:
            import nappy.pmd.pmd_wrapper  # noqa: F401
        except ImportError:
            self.skipTest('Fortran pmd_wrapper not compiled; skipping cross-validation')

        # scipy path
        orig = sys.modules.get('nappy.pmd.pairlist')
        sys.modules['nappy.pmd.pairlist'] = None
        try:
            nsys_s = _make_fcc(a=4.0, n=3)
            nsys_s.make_pair_list(rcut=3.2)
            nbs_s = [sorted(nsys_s.atoms['neighbors'][i])
                     for i in range(nsys_s.num_atoms())]
        finally:
            sys.modules['nappy.pmd.pairlist'] = orig

        # Fortran path
        nsys_f = _make_fcc(a=4.0, n=3)
        nsys_f.make_pair_list(rcut=3.2)
        nbs_f = [sorted(nsys_f.atoms['neighbors'][i])
                 for i in range(nsys_f.num_atoms())]

        for ia, (s, f) in enumerate(zip(nbs_s, nbs_f)):
            self.assertEqual(s, f, f'Atom {ia}: scipy={s[:3]}... fortran={f[:3]}...')


if __name__ == '__main__':
    unittest.main()
