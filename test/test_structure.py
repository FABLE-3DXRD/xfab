import unittest
from CifFile import ReadCif
import numpy as n

from xfab import structure

class test_cifread(unittest.TestCase):
    def test_cifopen(self):  ## test method names begin 'test*'
        mylist =  structure.build_atomlist()
        mylist.CIFopen('PPA.cif','oPPA')
    def test_cifread(self):
        mylist =  structure.build_atomlist()
        mylist.CIFread('PPA.cif','oPPA')
        self.assertEquals([8.5312,4.8321,10.125,90.00,92.031,90.00],
                          mylist.atomlist.cell)

class test_structurefactor(unittest.TestCase):
    def test_sfcalc(self):  ## test method names begin 'test*'
        # Read the cif
        mylist =  structure.build_atomlist()
        mylist.CIFread('PPA.cif','oPPA')

        # Read the fcf
        fcf = ReadCif('oPPA.fcf')['oPPA']
        for i in range(500):# Consider only the first 500 reflections
            hkl =[eval(fcf['_refln_index_h'][i]),
                  eval(fcf['_refln_index_k'][i]),
                  eval(fcf['_refln_index_l'][i])]
            (Fr, Fi) = structure.StructureFactor(hkl, mylist.atomlist.cell,
                                                 mylist.atomlist.sgname,
                                                 mylist.atomlist.atom,
                                                 mylist.atomlist.dispersion)
            F2 = Fr**2 + Fi**2
            reldif = F2/eval(fcf['_refln_F_squared_calc'][i])-1
            #print i, reldif, F2
            if F2 > 10: # Only compare those with an F^2 larger than ten 
                        # to avoid that the very weak reflections which
                        # have a relative difference that are slightly larger
                self.assertAlmostEquals(reldif,0,2)

if __name__ == '__main__':
    unittest.main()
