import unittest
import numpy as n
from xfab import tools


class test_euler2u(unittest.TestCase):
    def test1(self):  
        phi1 = 0.0
        PHI  = 0.0
        phi2 = 0.0
        Umat = tools.euler2U(phi1,PHI,phi2)
        diff = n.abs(Umat-n.eye(3)).sum()
        self.assertAlmostEquals(diff,0,9)

    def test2(self):  
        phi1 = 0.1
        PHI  = 0.0
        phi2 = -0.1
        Umat = tools.euler2U(phi1,PHI,phi2)
        diff = n.abs(Umat-n.eye(3)).sum()
        self.assertAlmostEquals(diff,0,9)

class test_rodrigues(unittest.TestCase):

    def test_rod2U2rod(self):  
        rodvec = n.array([0.23,-0.34,0.7])
        Umat = tools.rod2U(rodvec)
        rodvec2 = tools.U2rod(Umat)

        diff = n.abs(rodvec-rodvec2).sum()
        self.assertAlmostEquals(diff,0,9)

        
    def test_U2rod2U(self):  
        phi1 = 0.2
        PHI  = 0.4
        phi2 = 0.1
        Umat = tools.euler2U(phi1,PHI,phi2)
        rodvec = tools.U2rod(Umat)
        Umat2 = tools.rod2U(rodvec)

        diff = n.abs(Umat-Umat2).sum()
        self.assertAlmostEquals(diff,0,9)


    def test_ubi2rod(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        cell = [3,4,5,80,95,100]
        Umat = tools.euler2U(phi1,PHI,phi2)
        ubi = n.linalg.inv(n.dot(Umat, tools.FormB(cell)))*2*n.pi
        rodubi = tools.ubi2rod(ubi)
        rodU = tools.U2rod(Umat)
        diff = n.abs(rodubi-rodU).sum()
        print diff
        self.assertAlmostEquals(diff,0,9)


if __name__ == '__main__':
    unittest.main()
