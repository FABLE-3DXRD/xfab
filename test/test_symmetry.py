from __future__ import absolute_import
from __future__ import print_function
import unittest
import numpy as np

from xfab import symmetry
from six.moves import range


class test_Umis(unittest.TestCase):
    def setUp(self):
        np.random.seed(46)

    def test_identical_rotations(self):
        for i in range(1, 8):
            m = symmetry.Umis(np.eye(3, 3), np.eye(3, 3), crystal_system=i)
            self.assertAlmostEqual(np.min(m), 0)

    def test_small_perturbation_rotations(self):
        for _ in range(100):
            angle = 2 * (np.random.rand() - 0.5) * 180
            dangle = 2 * (np.random.rand() - 0.5)

            U1 = self._get_rot_z_mat(angle)
            U2 = self._get_rot_z_mat(angle + dangle)

            for i in range(1, 8):
                m = symmetry.Umis(U1, U2, crystal_system=i)
                self.assertLessEqual(np.abs(np.min(m[:, 1]) - np.abs(dangle)), 1e-6)

    def test_90dgr_flipped_cubic(self):
        for _ in range(100):
            angle = 2 * (np.random.rand() - 0.5) * 180
            dangle = 90

            U1 = self._get_rot_z_mat(angle)
            U2 = self._get_rot_z_mat(angle + dangle)

            m = symmetry.Umis(U1, U2, crystal_system=7)
            self.assertAlmostEqual(np.min(m[:, 1]), 0)
            self.assertLessEqual(np.min(np.abs(m[:, 1] - 90)), 1e-6)

    def test_old_Umis(self):
        for _ in range(100):
            U1 = self.random_uniform_orientation()
            U2 = self.random_uniform_orientation()
            for i in range(1, 8):
                m1 = self._Umis_old(U1, U2, crystal_system=i)
                m2 = symmetry.Umis(U1, U2, crystal_system=i)
                self.assertLessEqual(np.max(np.abs(m1 - m2)), 1e-6)

    def _get_rot_z_mat(self, angle):
        c, s = np.cos(np.radians(angle)), np.sin(np.radians(angle))
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def _Umis_old(self, umat_1, umat_2, crystal_system):
        """Old verified implementation (prior to 2022 implementation of Umis)"""
        rot = symmetry.rotations(crystal_system)
        t_save = np.zeros((0, 2))
        for k in range(rot.shape[0]):
            g_vector = np.dot(umat_2, np.transpose(np.dot(umat_1, rot[k])))
            length = (g_vector.diagonal().sum() - 1.0) / 2.0
            if abs(length) > 1.00000000:
                if length > 1:
                    length = 1.0
                else:
                    length = -1.0
            misangle = np.arccos(length)
            t_save = np.concatenate((t_save, np.array([[k, misangle * 180 / np.pi]])), 0)
        return t_save

    def random_uniform_orientation(self):
        """Generate random rotations using QR.
        see M. Ozols, “How to generate a random unitary matrix,” 2009.
        """
        Z = np.random.standard_normal((3, 3))
        Q, _ = np.linalg.qr(Z)
        assert np.sum(np.abs(np.eye(3, 3) - Q.T.dot(Q))) < 1e-8
        assert np.abs(np.linalg.det(Q) - 1.0) < 1e-8
        return Q


if __name__ == "__main__":
    unittest.main()
