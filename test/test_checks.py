import numpy as np
import xfab
from xfab import tools
import unittest


class test_checks(unittest.TestCase):
    def setUp(self):
        np.random.seed(87)

    def test_default_activated(self):
        tools.u_to_euler(np.random.rand(3, 3))

    def test_deactivate(self):
        xfab.CHECKS.activated = False
        tools.u_to_euler(np.random.rand(3, 3))

    def test_reactivate(self):
        xfab.CHECKS.activated = True
        try:
            tools.u_to_euler(np.random.rand(3, 3))
        except ValueError as e:
            self.assertEqual(str(e), "orientation matrix U is not unitary, np.dot(U.T, U)!=np.eye(3,3)")


if __name__ == "__main__":
    unittest.main()
