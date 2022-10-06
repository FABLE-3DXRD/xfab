from __future__ import absolute_import
import unittest
import numpy as np
from xfab import laue
from six.moves import range

np.random.seed(0) # to make all unittest repeatable


class test_euler2u(unittest.TestCase):
    def test1(self):  
        phi1 = 0.0
        PHI  = 0.0
        phi2 = 0.0
        Umat = laue.euler_to_u(phi1,PHI,phi2)
        diff = np.abs(Umat-np.eye(3)).sum()
        self.assertAlmostEqual(diff,0,9)

    def test2(self):  
        phi1 = 0.1
        PHI  = 0.0
        phi2 = -0.1 + 2*np.pi
        Umat = laue.euler_to_u(phi1,PHI,phi2)
        diff = np.abs(Umat-np.eye(3)).sum()
        self.assertAlmostEqual(diff,0,9)

    def test_gridded_euler_space(self):
        '''Test that euler_to_u() and u_to_euler() give consistent output
        for a coarsly gridded Euler space, phi1,PHI,phi2 = [0 , 2*pi].
        between the two functions.
        '''
        gridded_angles = np.linspace(0, 2*np.pi, 20)
        self._check_angles( gridded_angles )

    def test_special_euler_angles(self):
        '''Test that euler_to_u() and u_to_euler() give consistent output
        for special angles such as: pi, pi/2, 0 ... etc
        '''
        special_angles = np.array([0, np.pi/2, np.pi, 2*np.pi])
        self._check_angles( special_angles )

    def _check_angles(self, angles):
        '''Check consistency for euler_to_u() and u_to_euler(). I.e check that 
        the orientation matrix constructed for some set of angles is recovered
        again when passing between the two functions.
        '''
        for phi1 in angles:
            for PHI in angles:
                for phi2 in angles:
                    U = laue.euler_to_u(phi1, PHI, phi2)
                    phi1_new, PHI_new, phi2_new = laue.u_to_euler( U )
                    U_new = laue.euler_to_u(phi1_new, PHI_new, phi2_new)
                    maxdiff = np.max( np.abs( U - U_new ) )
                    self.assertTrue( maxdiff < 1e-8 )
                    self.assertTrue( 0 <= phi1_new <= 2*np.pi )
                    self.assertTrue( 0 <= PHI_new  <= 2*np.pi )
                    self.assertTrue( 0 <= phi2_new <= 2*np.pi )

class test_rodrigues(unittest.TestCase):

    def test_rod2U2rod(self):  
        rodvec = np.array([0.23,-0.34,0.7])
        Umat = laue.rod_to_u(rodvec)
        rodvec2 = laue.u_to_rod(Umat)

        diff = np.abs(rodvec-rodvec2).sum()
        self.assertAlmostEqual(diff,0,9)

        
    def test_U2rod2U(self):  
        phi1 = 0.2
        PHI  = 0.4
        phi2 = 0.1
        Umat = laue.euler_to_u(phi1,PHI,phi2)
        rodvec = laue.u_to_rod(Umat)
        Umat2 = laue.rod_to_u(rodvec)

        diff = np.abs(Umat-Umat2).sum()
        self.assertAlmostEqual(diff,0,9)


    def test_ubi2rod(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        cell = [3,4,5,80,95,100]
        Umat = laue.euler_to_u(phi1,PHI,phi2)
        ubi = np.linalg.inv(np.dot(Umat, laue.form_b_mat(cell)))
        rodubi = laue.ubi_to_rod(ubi)
        rodU = laue.u_to_rod(Umat)
        diff = np.abs(rodubi-rodU).sum()
        self.assertAlmostEqual(diff,0,9)

class test_u2ubi(unittest.TestCase):
    def test1(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        cell = [3,4,5,80,95,100]
        Umat = laue.euler_to_u(phi1,PHI,phi2)
        Bmat = laue.form_b_mat(cell)
        ubi = np.linalg.inv(np.dot(Umat, Bmat))
        ubi2 = laue.u_to_ubi(Umat,cell)
        (U2,B2) = laue.ubi_to_u_b(ubi2)
        diff = np.abs(ubi-ubi2).sum()
        self.assertAlmostEqual(diff,0,9)  
        diff = np.abs(Umat-U2).sum()
        self.assertAlmostEqual(diff,0,9)  
        diff = np.abs(Bmat-B2).sum()
        self.assertAlmostEqual(diff,0,9)  

    def test_precision_lost_in_ubi_and_b_transforms(self):
        # Test that small errors are allowed by xfab.checks regardless of dtypes.
        unit_cell = [3, 4, 5, 80, 95, 100]
        Bmat = laue.form_b_mat(unit_cell)
        for _ in range(10):
            Umat, _ = np.linalg.qr( np.random.standard_normal((3, 3)) )
            ubi = np.linalg.inv(np.dot(Umat, Bmat))
            for dtype in [float, np.float64, np.float32, np.float16]:
                _ = laue.ubi_to_u(ubi.copy().astype(dtype))
                _, _ = laue.ubi_to_u_and_eps(ubi.copy().astype(dtype), unit_cell)
                _ = laue.b_to_cell(Bmat.copy().astype(dtype))
                _ = laue.b_to_epsilon_old(Bmat.copy().astype(dtype), unit_cell)
                _ = laue.b_to_epsilon(Bmat.copy().astype(dtype), unit_cell)

class test_ubi_to_u_and_eps(unittest.TestCase):
    def test1(self):
        unit_cell = [1., 1., 1., 90., 90., 90.]
        eps_true = [0.01, 0.0, 0.024, -0.03, 0.3, 0.0]
        b_true = laue.epsilon_to_b(eps_true, unit_cell)
        u_true = np.eye(3)
        ubi = np.linalg.inv( u_true.dot(b_true) )
        u, eps = laue.ubi_to_u_and_eps(ubi, unit_cell)
        for e1,e2 in zip(eps,eps_true):
            self.assertAlmostEqual(e1,e2)
        for u1,u2 in zip(u_true.flatten(), u.flatten()):
            self.assertAlmostEqual(u1,u2)

class test_twotheta(unittest.TestCase):

    def test_tth(self):
        # generate random gvector
        hkl = np.array([round(np.random.rand()*10-5),round(np.random.rand()*10-5),round(np.random.rand()*10-5)])
        ucell = np.array([3.5+np.random.rand(),3.5+np.random.rand(),3.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand()])
        B = laue.form_b_mat(ucell)
        U = laue.euler_to_u(np.random.rand()*2.*np.pi,np.random.rand()*2.*np.pi,np.random.rand()*np.pi)
        wavelength = 0.95 + np.random.rand()*0.1
        gvec = np.dot(U,np.dot(B,hkl))
        tth = laue.tth(ucell, hkl, wavelength)
        tth2 = laue.tth2(gvec,wavelength)
        diff = np.abs(tth-tth2)
        self.assertAlmostEqual(diff,0,9)

class test_ubi_to_cell(unittest.TestCase):

    def test1(self):
        unit_cell_1  =  [2., 3., 4., 90.,90.,120.]
        B = laue.form_b_mat(unit_cell_1)
        U = np.eye(3, 3, dtype=np.float64)
        ubi = np.linalg.inv( U.dot(B) )
        unit_cell_2 = laue.ubi_to_cell(ubi)
        for c1,c2 in zip(unit_cell_1, unit_cell_2):
            self.assertAlmostEqual( c1, c2, msg='unit cell is not preserved over B matrix cycling')

class test_general_orientation(unittest.TestCase):

    def test_find_omega_general_nowedge(self):
        # generate random gvector
        hkl = np.array([round(np.random.rand()*10-5),round(np.random.rand()*10-5),round(np.random.rand()*10-5)])
        ucell = np.array([3.5+np.random.rand(),3.5+np.random.rand(),3.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand()])
        B = laue.form_b_mat(ucell)
        U = laue.euler_to_u(np.random.rand()*2.*np.pi,np.random.rand()*2.*np.pi,np.random.rand()*np.pi)
        wavelength = 0.95 + np.random.rand()*0.1
        gvec = np.dot(U,np.dot(B,hkl))
        tth = laue.tth(ucell, hkl, wavelength)
        # calculate corresponding eta and Omega using laue.find_omega_general
        (omega1, eta1) = laue.find_omega_general(gvec*wavelength/(4.*np.pi),tth,0,0)
        Om1 = []
        for i in range(len(omega1)):
            Om1.append(laue.form_omega_mat_general(omega1[i],0,0))  
        # calculate corresponding eta and Omega using laue.find_omega_wedge
        omega2 = laue.find_omega(gvec,tth)
        Om2 = []
        for i in range(len(omega2)):
            Om2.append(laue.form_omega_mat(omega2[i]))
        #assert  
        for i in range(len(omega1)):    
#            print Om1[i]
#            print Om2[i]
            diff = np.abs(Om1[i]-Om2[i]).sum()
            self.assertAlmostEqual(diff,0,9)        



class test_ABepsilon(unittest.TestCase):

    def test_ucell2A2ucell(self):
        ucell = np.array([3.5+np.random.rand(),3.5+np.random.rand(),3.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand()])
        A = laue.form_a_mat(ucell)
        ucell2 = laue.a_to_cell(A)
        diff = np.abs(ucell-ucell2).sum()
        self.assertAlmostEqual(diff,0,9)        
        
        
    def test_ucell2B2ucell(self):
        ucell = np.array([3.5+np.random.rand(),3.5+np.random.rand(),3.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand()])
        B = laue.form_b_mat(ucell)
        ucell2 = laue.b_to_cell(B)
        diff = np.abs(ucell-ucell2).sum()
        self.assertAlmostEqual(diff,0,9)  

    def test_epsilon2B2epsilon(self):
        ucell = np.array([3.5+np.random.rand(),3.5+np.random.rand(),3.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand()])
        eps = (np.random.rand(6)-.5)/1000.
        B = laue.epsilon_to_b(eps,ucell)
        eps2 = laue.b_to_epsilon(B,ucell)
        
class test_qr(unittest.TestCase):

    def test_UdotBtoUandB(self):
        ucell = np.array([3.5+np.random.rand(),3.5+np.random.rand(),3.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand(),89.5+np.random.rand()])
        eps = (np.random.rand(6)-.5)/1000.
        B = laue.epsilon_to_b(eps,ucell)
        U = laue.euler_to_u(np.random.rand()*2.*np.pi,np.random.rand()*2.*np.pi,np.random.rand()*np.pi)
        UB = np.dot(U,B)
        (U1,B1) = laue.ub_to_u_b(UB)
        diffU = np.abs(U-U1).sum()
        diffB = np.abs(B-B1).sum()
        self.assertAlmostEqual(diffU,0,9)  
        self.assertAlmostEqual(diffB,0,9)  
        

    def test_ubi2rod(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        cell = [3,4,5,80,95,100]
        Umat = laue.euler_to_u(phi1,PHI,phi2)
        Bmat = laue.form_b_mat(cell)
        ubi = np.linalg.inv(np.dot(Umat, Bmat))
        (U1,B1) = laue.ubi_to_u_b(ubi)
        diffU = np.abs(Umat-U1).sum()
        diffB = np.abs(Bmat-B1).sum()
        self.assertAlmostEqual(diffU,0,9)  
        self.assertAlmostEqual(diffB,0,9)  

class test_reduce_cell(unittest.TestCase):
    def test_1(self):
        cell = np.array([9.07599708738,
                        6.05007626616,
                        44.33571511, 
                        97.838350766558762, 
                        90.0, 
                        90.0])

        red_cell = laue.reduce_cell(cell)
        
        red_cell_ref = np.array([9.0759970873799993, 
                                6.0500773546922257,
                                43.921476668199631,
                                89.965630135621325,
                                89.999999999999986,
                                90.0])

        diff_cell = np.abs(red_cell - red_cell_ref)
        for diff_par in diff_cell:
            self.assertAlmostEqual(diff_par,0,9)  
    
if __name__ == '__main__':
    unittest.main()
