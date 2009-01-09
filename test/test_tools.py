import unittest
import numpy as n
from xfab import tools


class test_euler2u(unittest.TestCase):
    def test1(self):  
        phi1 = 0.0
        PHI  = 0.0
        phi2 = 0.0
        Umat = tools.euler2U(phi1,PHI,phi2)
        self.assertEquals(Umat.ravel().all(),n.eye(3).ravel().all())

    def test2(self):  
        phi1 = 0.1
        PHI  = 0.0
        phi2 = -0.1
        Umat = tools.euler2U(phi1,PHI,phi2)
        self.assertEquals(Umat.ravel().all(),n.eye(3).ravel().all())

class test_rodrigues(unittest.TestCase):
    def test_compare_U2rod_methods(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        Umat = tools.euler2U(phi1,PHI,phi2)
        rodvec = tools.U2rod(Umat)
        rodvec_old = tools.U2rod_old(Umat)
        self.assertEquals(rodvec.ravel().all(),rodvec_old.ravel().all())

    def test_rod2U2rod(self):  
        rodvec = n.array([0.23,-0.34,0.7])
        Umat = tools.rod2U(rodvec)
        rodvec2 = tools.U2rod(Umat)
        self.assertEquals(rodvec.ravel().all(),rodvec2.ravel().all())
        
    def test_U2rod2U(self):  
        phi1 = 0.2
        PHI  = 0.4
        phi2 = 0.1
        Umat = tools.euler2U(phi1,PHI,phi2)
        rodvec = tools.U2rod(Umat)
        Umat2 = tools.rod2U(rodvec)
        self.assertEquals(Umat.ravel().all(),Umat2.ravel().all())

if __name__ == '__main__':
    unittest.main()
