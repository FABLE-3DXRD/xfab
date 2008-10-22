import unittest
import numpy as n
from xfab import tools,detector

class test_detector_flips(unittest.TestCase):

    def test1(self):  ## o11, o12, o21, o22 = 1, 0, 0, 1
        a = n.arange(25).reshape(5,5)
        b = n.transpose(a)
        b2 = detector.trans_orientation(a,1,0,0,1)
        self.assertEquals(b.tolist(),b2.tolist())

    def test1_inverse(self):   ## o11, o12, o21, o22 = 1, 0, 0, 1
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,1,0,0,1)
        b = detector.trans_orientation(b2,1,0,0,1,'inverse')
        self.assertEquals(a.tolist(),b.tolist())

    def test2(self):  ## o11, o12, o21, o22 = -1, 0, 0, 1
        a = n.arange(25).reshape(5,5)
        b = n.fliplr(n.transpose(a))
        b2 = detector.trans_orientation(a,-1,0,0,1)
        self.assertEquals(b.tolist(),b2.tolist())

    def test2_inverse(self):   ## o11, o12, o21, o22 = -1, 0, 0, 1
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,-1,0,0,1)
        b = detector.trans_orientation(b2,-1,0,0,1,'inverse')
        self.assertEquals(a.tolist(),b.tolist())

    def test3(self):  ## test flip matching frelon2k/4m
        ## o11, o12, o21, o22 = 1, 0, 0, -1
        a = n.arange(25).reshape(5,5)
        b = n.flipud(n.transpose(a))
        b2 = detector.trans_orientation(a,1,0,0,-1)
        self.assertEquals(b.tolist(),b2.tolist())

    def test3_inverse(self):  ## test flip matching frelon2k/4m
        ## o11, o12, o21, o22 = 1, 0, 0, -1
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,1,0,0,-1)
        b = detector.trans_orientation(b2,1,0,0,-1,'inverse')
        self.assertEquals(a.tolist(),b.tolist())

    def test4(self): ## o11, o12, o21, o22 = -1, 0, 0, -1
        a = n.arange(25).reshape(5,5)
        b = n.fliplr(n.flipud(n.transpose(a)))
        b2 = detector.trans_orientation(a,-1,0,0,-1)
        self.assertEquals(b.tolist(),b2.tolist())

    def test4_inverse(self): ## o11, o12, o21, o22 = -1, 0, 0, -1
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,-1,0,0,-1)
        b = detector.trans_orientation(b2,-1,0,0,-1,'inverse')
        self.assertEquals(a.tolist(),b.tolist())

    def test5(self): ## o11, o12, o21, o22 = 0, 1, 1, 0
        a = n.arange(25).reshape(5,5)
        b = a
        b2 = detector.trans_orientation(a,0,1,1,0)
        self.assertEquals(b.tolist(),b2.tolist())

    def test5_inverse(self):  ## o11, o12, o21, o22 = 0, 1, 1, 0
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,0,1,1,0)
        b = detector.trans_orientation(b2,0,1,1,0,'inverse')
        self.assertEquals(a.tolist(),b.tolist())

    def test6(self):  ## o11, o12, o21, o22 = 0, -1, 1, 0
        a = n.arange(25).reshape(5,5)
        b = n.fliplr(a)
        b2 = detector.trans_orientation(a,0,-1,1,0)
        self.assertEquals(b.tolist(),b2.tolist())

    def test6_inverse(self):   ## o11, o12, o21, o22 = 0, -1, 1, 0
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,0,-1,1,0)
        b = detector.trans_orientation(b2,0,-1,1,0,'inverse')
        self.assertEquals(a.tolist(),b.tolist())

    def test7(self):   ## o11, o12, o21, o22 = 0, 1, -1, 0
        a = n.arange(25).reshape(5,5)
        b = n.flipud(a)
        b2 = detector.trans_orientation(a,0,1,-1,0)
        self.assertEquals(b.tolist(),b2.tolist())

    def test7_inverse(self):     ## o11, o12, o21, o22 = 0, 1, -1, 0
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,0,1,-1,0)
        b = detector.trans_orientation(b2,0,1,-1,0,'inverse')
        self.assertEquals(a.tolist(),b.tolist())

    def test8(self):     ## o11, o12, o21, o22 = 0, -1, -1, 0 
        a = n.arange(25).reshape(5,5)
        b = n.fliplr(n.flipud(a))
        b2 = detector.trans_orientation(a,0,-1,-1,0)
        self.assertEquals(b.tolist(),b2.tolist())

    def test8_inverse(self):     ## o11, o12, o21, o22 = 0, -1, -1, 0
        a = n.arange(25).reshape(5,5)
        b2 = detector.trans_orientation(a,0,-1,-1,0)
        b = detector.trans_orientation(b2,0,-1,-1,0,'inverse')
        self.assertEquals(a.tolist(),b.tolist())


class test_detector_coord_transform(unittest.TestCase):
    def test1(self):  ## o11, o12, o21, o22 = 1, 0, 0, 1
        xy = [10,20]
        (dety,detz) = detector.detyz2xy(xy,1,0,0,1,1024,1024)
        self.assertEquals(xy,[detz,dety])

    def test2(self):  ## o11, o12, o21, o22 = 1, 0, 0, 1
        detyz = [10,20]
        (x,y) = detector.xy2detyz(detyz,1,0,0,1,1024,1024)
        self.assertEquals(detyz,[y,x])

    def test3(self):  ## o11, o12, o21, o22 = -1, 0, 0, -1
        detyz = n.array([841.38745747, 62.2754412563])
        xy = detector.xy2detyz(detyz,-1,0,0,-1,1024,1024)
        detyz2 = detector.detyz2xy(xy,-1,0,0,-1,1024,1024)
        self.assertEquals(detyz.all(),detyz2.all())

class test_detector_coord(unittest.TestCase):
    def test_compare_calc(self):  
        tth =  21.836334 *n.pi/180.
        omega = -89.247851 *n.pi/180.
        eta =  45.480703*n.pi/180.
        (dety_orig, detz_orig) = (109.418256, 925.772491)
        gv = n.dot(tools.OMEGA(omega),n.array([3.260078, -0.928075, 3.217533]))
        wavelength = 0.5092836 # in angstrom
        distance  = 135.00                        # sample-detector distance (mm)
        dety_center = 521.5                              # beamcenter, y in pixel coordinatees
        detz_center = 531.5                              # beamcenter, z in pixel coordinatees
        y_size = 0.0936   # Pixel size y (mm)
        z_size = 0.0962   # Pixel size z (mm)
        dety_size = 1024.0     # detector y size (pixels)
        detz_size = 1024.0     # detector z size (pixels)
        tilt_x =   0.0        # detector tilt counterclockwise around lab x axis in rad 
        tilt_y =   0.0        # detector tilt counterclockwise around lab y axis in rad 
        tilt_z =   0.0       # detector tilt counterclockwise around lab z axis in rad 
        tx = 0
        ty = 0
        tz = 0
        R_tilt = tools.detect_tilt(tilt_x,tilt_y,tilt_z)
        costth = n.cos(tth)

        (dety1,detz1) = detector.det_coor(gv,costth , wavelength, distance, y_size, z_size, dety_center, detz_center,
             R_tilt, tx, ty, tz)
        #print dety1,detz1
        (dety2,detz2) = detector.det_coor2(tth, eta, distance, y_size, z_size, dety_center, detz_center,
             R_tilt, tx, ty, tz)
        #print dety2,detz2

        self.assertAlmostEquals(dety1,dety2,4)
        self.assertAlmostEquals(detz1,detz2,4)

if __name__ == '__main__':
    unittest.main()