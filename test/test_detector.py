import unittest
import numpy as n
from xfab import detector

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


if __name__ == '__main__':
    unittest.main()
