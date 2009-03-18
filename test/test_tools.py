import unittest
import numpy as n
from xfab import tools


class test_euler2u(unittest.TestCase):
    def test1(self):  
        phi1 = 0.0
        PHI  = 0.0
        phi2 = 0.0
        Umat = tools.euler_to_u(phi1,PHI,phi2)
        diff = n.abs(Umat-n.eye(3)).sum()
        self.assertAlmostEquals(diff,0,9)

    def test2(self):  
        phi1 = 0.1
        PHI  = 0.0
        phi2 = -0.1
        Umat = tools.euler_to_u(phi1,PHI,phi2)
        diff = n.abs(Umat-n.eye(3)).sum()
        self.assertAlmostEquals(diff,0,9)

class test_rodrigues(unittest.TestCase):

    def test_rod2U2rod(self):  
        rodvec = n.array([0.23,-0.34,0.7])
        Umat = tools.rod_to_u(rodvec)
        rodvec2 = tools.u_to_rod(Umat)

        diff = n.abs(rodvec-rodvec2).sum()
        self.assertAlmostEquals(diff,0,9)

        
    def test_U2rod2U(self):  
        phi1 = 0.2
        PHI  = 0.4
        phi2 = 0.1
        Umat = tools.euler_to_u(phi1,PHI,phi2)
        rodvec = tools.u_to_rod(Umat)
        Umat2 = tools.rod_to_u(rodvec)

        diff = n.abs(Umat-Umat2).sum()
        self.assertAlmostEquals(diff,0,9)


    def test_ubi2rod(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        cell = [3,4,5,80,95,100]
        Umat = tools.euler_to_u(phi1,PHI,phi2)
        ubi = n.linalg.inv(n.dot(Umat, tools.form_b_mat(cell)))*2*n.pi
        rodubi = tools.ubi_to_rod(ubi)
        rodU = tools.u_to_rod(Umat)
        diff = n.abs(rodubi-rodU).sum()
        self.assertAlmostEquals(diff,0,9)

class test_u2ubi(unittest.TestCase):
    def test1(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        cell = [3,4,5,80,95,100]
        Umat = tools.euler_to_u(phi1,PHI,phi2)
        Bmat = tools.form_b_mat(cell)
        ubi = n.linalg.inv(n.dot(Umat, Bmat))*2*n.pi
        ubi2 = tools.u_to_ubi(Umat,cell)
        (U2,B2) = tools.ubi_to_u_b(ubi2)
        diff = n.abs(ubi-ubi2).sum()
        self.assertAlmostEquals(diff,0,9)  
        diff = n.abs(Umat-U2).sum()
        self.assertAlmostEquals(diff,0,9)  
        diff = n.abs(Bmat-B2).sum()
        self.assertAlmostEquals(diff,0,9)  



class test_twotheta(unittest.TestCase):

    def test_tth(self):
        # generate random gvector
        hkl = n.array([round(n.random.rand()*10-5),round(n.random.rand()*10-5),round(n.random.rand()*10-5)])
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        B = tools.form_b_mat(ucell)
        U = tools.euler_to_u(n.random.rand()*2.*n.pi,n.random.rand()*2.*n.pi,n.random.rand()*n.pi)
        wavelength = 0.95 + n.random.rand()*0.1
        gvec = n.dot(U,n.dot(B,hkl))
        tth = tools.tth(ucell, hkl, wavelength)
        tth2 = tools.tth2(gvec,wavelength)
        diff = n.abs(tth-tth2)
        self.assertAlmostEquals(diff,0,9)
        
class test_quarternions(unittest.TestCase):

    def test_quart2Omega(self):
        omega = n.random.rand()*360.
        Omega = tools.form_omega_mat(omega*n.pi/180.)
        Omega_quart = tools.quart_to_omega(omega,0,0)
        diff = n.abs(Omega-Omega_quart).sum()
        self.assertAlmostEquals(diff,0,9)
        
    def test_find_omega(self):
        # generate random gvector
        hkl = n.array([round(n.random.rand()*10-5),round(n.random.rand()*10-5),round(n.random.rand()*10-5)])
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        B = tools.form_b_mat(ucell)
        U = tools.euler_to_u(n.random.rand()*2.*n.pi,n.random.rand()*2.*n.pi,n.random.rand()*n.pi)
        wavelength = 0.95 + n.random.rand()*0.1
        gvec = n.dot(U,n.dot(B,hkl))
        tth = tools.tth(ucell, hkl, wavelength)        
        # calculate corresponding eta and Omega using tools.find_omega_quart
        (omega1, eta1) = tools.find_omega_quart(gvec*wavelength/(4.*n.pi),tth,0,0)
        Om1 = []
        for i in range(len(omega1)):
            Om1.append(tools.quart_to_omega(omega1[i]*180./n.pi,0,0))  
        # calculate corresponding eta and Omega using tools.find_omega
        omega2 = tools.find_omega(gvec,tth)
        Om2 = []
        for i in range(len(omega2)):
            Om2.append(tools.form_omega_mat(omega2[i]))
        # assert  
        for i in range(len(eta1)):          
            diff = n.abs(Om1[i]-Om2[i]).sum()
            self.assertAlmostEquals(diff,0,9)        
            
    def test_find_omega_wedge(self):
        # generate random gvector
        hkl = n.array([round(n.random.rand()*10-5),round(n.random.rand()*10-5),round(n.random.rand()*10-5)])
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        B = tools.form_b_mat(ucell)
        U = tools.euler_to_u(n.random.rand()*2.*n.pi,n.random.rand()*2.*n.pi,n.random.rand()*n.pi)
        wavelength = 0.95 + n.random.rand()*0.1
        gvec = n.dot(U,n.dot(B,hkl))
        tth = tools.tth(ucell, hkl, wavelength)
        wedge = n.random.rand()-.5
        # calculate corresponding eta and Omega using tools.find_omega_quart
        (omega1, eta1) = tools.find_omega_quart(gvec*wavelength/(4.*n.pi),tth,0,wedge)
        Om1 = []
        for i in range(len(omega1)):
            Om1.append(tools.quart_to_omega(omega1[i]*180./n.pi,0,wedge))  
        # calculate corresponding eta and Omega using tools.find_omega_wedge
        Phi_y = n.array([[ n.cos(wedge), 0, n.sin(wedge)],
                         [0         , 1, 0        ],
                         [-n.sin(wedge), 0, n.cos(wedge)]])
        gvec = n.dot(n.transpose(Phi_y),gvec)
        (omega2, eta2) = tools.find_omega_wedge(gvec,tth,-wedge)
        Om2 = []
        for i in range(len(omega2)):
            Om2.append(n.dot(Phi_y,n.dot(tools.form_omega_mat(omega2[i]),n.transpose(Phi_y))))
        #assert  
        for i in range(len(eta1)):          
            diff = n.abs(eta1[i]-eta2[i])
            self.assertAlmostEquals(diff,0,9)
            diff = n.abs(Om1[i]-Om2[i]).sum()
            self.assertAlmostEquals(diff,0,9)        

class test_general_orientation(unittest.TestCase):

    def test_find_omega_general_nowedge(self):
        # generate random gvector
        hkl = n.array([round(n.random.rand()*10-5),round(n.random.rand()*10-5),round(n.random.rand()*10-5)])
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        B = tools.form_b_mat(ucell)
        U = tools.euler_to_u(n.random.rand()*2.*n.pi,n.random.rand()*2.*n.pi,n.random.rand()*n.pi)
        wavelength = 0.95 + n.random.rand()*0.1
        gvec = n.dot(U,n.dot(B,hkl))
        tth = tools.tth(ucell, hkl, wavelength)
        # calculate corresponding eta and Omega using tools.find_omega_general
        (omega1, eta1) = tools.find_omega_general(gvec*wavelength/(4.*n.pi),tth,0,0)
        Om1 = []
        for i in range(len(omega1)):
            Om1.append(tools.form_omega_mat_general(omega1[i],0,0))  
        # calculate corresponding eta and Omega using tools.find_omega_wedge
        omega2 = tools.find_omega(gvec,tth)
        Om2 = []
        for i in range(len(omega2)):
            Om2.append(tools.form_omega_mat(omega2[i]))
        #assert  
        for i in range(len(omega1)):    
#            print Om1[i]
#            print Om2[i]
            diff = n.abs(Om1[i]-Om2[i]).sum()
            self.assertAlmostEquals(diff,0,9)        

    def test_find_omega_general_wedge(self):
        # generate random gvector
        hkl = n.array([round(n.random.rand()*10-5),round(n.random.rand()*10-5),round(n.random.rand()*10-5)])
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        B = tools.form_b_mat(ucell)
        U = tools.euler_to_u(n.random.rand()*2.*n.pi,n.random.rand()*2.*n.pi,n.random.rand()*n.pi)
        wavelength = 0.95 + n.random.rand()*0.1
        gvec = n.dot(U,n.dot(B,hkl))
        tth = tools.tth(ucell, hkl, wavelength)
#        wedge = n.random.rand()-.5
        wedge = n.random.rand()*10-5
        Phi_y = n.array([[ n.cos(wedge), 0, n.sin(wedge)],
                         [0         , 1, 0        ],
                         [-n.sin(wedge), 0, n.cos(wedge)]])
        # calculate corresponding eta and Omega using tools.find_omega_general
        (omega1, eta1) = tools.find_omega_general(gvec*wavelength/(4.*n.pi),tth,0,wedge)
        Om1 = []
        for i in range(len(omega1)):
            Om1.append(tools.form_omega_mat_general(omega1[i],0,wedge))  
        if len(eta1) == 2: 
            eta1 = n.array(eta1)
            Om1 = n.array(Om1)
            order = n.argsort(eta1)
            Om1 = Om1[order]
            eta1 = eta1[order]
        # calculate corresponding eta and Omega using tools.find_omega_wedge
        (omega2, eta2) = tools.find_omega_wedge(gvec,tth,-wedge)
        Om2 = []            
        for i in range(len(omega2)):
            Om2.append(n.dot(Phi_y,tools.form_omega_mat(omega2[i])))
        if len(eta2) == 2: 
            eta2 = n.array(eta2)
            Om2 = n.array(Om2)
            order = n.argsort(eta2)
            eta2 = eta2[order]
            Om2 = Om2[order]
        #assert  
        for i in range(len(omega1)):    
            if i==0:
                print 'test_find_omega_general_wedge'
#            print eta1[i],eta2[i]
#            print Om1[i]
#            print Om2[i]
            diff = n.abs(eta1[i]-eta2[i])
            self.assertAlmostEquals(diff,0,9)
            diff = n.abs(Om1[i]-Om2[i]).sum()
            self.assertAlmostEquals(diff,0,9)        


class test_ABepsilon(unittest.TestCase):

    def test_ucell2A2ucell(self):
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        A = tools.form_a_mat(ucell)
        ucell2 = tools.a_to_cell(A)
        diff = n.abs(ucell-ucell2).sum()
        self.assertAlmostEquals(diff,0,9)        
        
        
    def test_ucell2B2ucell(self):
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        B = tools.form_b_mat(ucell)
        ucell2 = tools.b_to_cell(B)
        diff = n.abs(ucell-ucell2).sum()
        self.assertAlmostEquals(diff,0,9)  

    def test_epsilon2B2epsilon(self):
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        eps = (n.random.rand(6)-.5)/1000.
        B = tools.epsilon_to_b(eps,ucell)
        eps2 = tools.b_to_epsilon(B,ucell)
        
class test_qr(unittest.TestCase):

    def test_UdotBtoUandB(self):
        ucell = n.array([3.5+n.random.rand(),3.5+n.random.rand(),3.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand(),89.5+n.random.rand()])
        eps = (n.random.rand(6)-.5)/1000.
        B = tools.epsilon_to_b(eps,ucell)
        U = tools.euler_to_u(n.random.rand()*2.*n.pi,n.random.rand()*2.*n.pi,n.random.rand()*n.pi)
        UB = n.dot(U,B)
        (U1,B1) = tools.ub_to_u_b(UB)
        diffU = n.abs(U-U1).sum()
        diffB = n.abs(B-B1).sum()
        self.assertAlmostEquals(diffU,0,9)  
        self.assertAlmostEquals(diffB,0,9)  
        

    def test_ubi2rod(self):  
        phi1 = 0.13
        PHI  = 0.4
        phi2 = 0.21
        cell = [3,4,5,80,95,100]
        Umat = tools.euler_to_u(phi1,PHI,phi2)
        Bmat = tools.form_b_mat(cell)
        ubi = n.linalg.inv(n.dot(Umat, Bmat))*2*n.pi
        (U1,B1) = tools.ubi_to_u_b(ubi)
        diffU = n.abs(Umat-U1).sum()
        diffB = n.abs(Bmat-B1).sum()
        self.assertAlmostEquals(diffU,0,9)  
        self.assertAlmostEquals(diffB,0,9)  

class test_reduce_cell(unittest.TestCase):
    def test_1(self):
        cell = n.array([9.07599708738,
                        6.05007626616,
                        44.33571511, 
                        97.838350766558762, 
                        90.0, 
                        90.0])

        red_cell = tools.reduce_cell(cell)
        
        red_cell_ref = n.array([9.0759970873799993, 
                                6.0500773546922257,
                                43.921476668199631,
                                89.965630135621325,
                                89.999999999999986,
                                90.0])

        diff_cell = n.abs(red_cell - red_cell_ref)
        for diff_par in diff_cell:
            self.assertAlmostEquals(diff_par,0,9)  
    
if __name__ == '__main__':
    unittest.main()
