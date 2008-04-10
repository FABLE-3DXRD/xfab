import unittest
from xfab import structure
class test_cifread(unittest.TestCase):
    def test_cifopen(self):  ## test method names begin 'test*'
        mylist =  structure.build_atomlist()
        mylist.CIFopen('PPA.cif','oPPA')
    def test_cifread(self):
        mylist =  structure.build_atomlist()
        mylist.CIFread('PPA.cif','oPPA')
        self.assertEquals([8.5312,4.8321,10.125,90.00,92.031,90.00], mylist.atomlist.cell)

if __name__ == '__main__':
    unittest.main()
