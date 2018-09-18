#!/usr/bin/env python3

import unittest
import sys
sys.path.insert(0,'..')
#from deviaTE import deviaTE_pileup as pileup
import deviaTE_pileup as pileup


class Test_is_snp(unittest.TestCase):
    
    def setUp(self):
        # set up site instance
        self.case1 = pileup.Site(pos=0, refbase='A', sid='a', fam='b')
        
    def test_is_snp1(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=0,
                                           C=0,
                                           G=0,
                                           T=0,
                                           cov=0), (False, False) )
        
    def test_is_snp2(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=100,
                                           C=0,
                                           G=0,
                                           T=0,
                                           cov=100), (False, False) )
        
    def test_is_snp3(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=0,
                                           C=100,
                                           G=0,
                                           T=0,
                                           cov=100), (False, True) )
        
    def test_is_snp4(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=90,
                                           C=10,
                                           G=0,
                                           T=0,
                                           cov=100), (True, False) )
        
    def test_is_snp5(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=10,
                                           C=90,
                                           G=0,
                                           T=0,
                                           cov=100), (True, False) )
        
    def test_is_snp6(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=50,
                                           C=50,
                                           G=0,
                                           T=0,
                                           cov=100), (True, False) )
        
    def test_is_snp7(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=33,
                                           C=33,
                                           G=33,
                                           T=0,
                                           cov=99), (True, False) )

    def test_is_snp8(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=25,
                                           C=25,
                                           G=25,
                                           T=25,
                                           cov=100), (True, False) )
        
    def test_is_snp9(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=1,
                                           C=99,
                                           G=0,
                                           T=0,
                                           cov=100), (True, False) )
        
    def test_is_snp10(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=99,
                                           C=1,
                                           G=0,
                                           T=0,
                                           cov=100), (False, False) )
        
    def test_is_snp11(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=9,
                                           C=1,
                                           G=0,
                                           T=0,
                                           cov=10), (False, False) )

    def test_is_snp12(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=1,
                                           C=9,
                                           G=0,
                                           T=0,
                                           cov=10), (True, False) )
        
    def test_is_snp13(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=99,
                                           C=1,
                                           G=0,
                                           T=0,
                                           cov=100), (False, False) )
        
    def test_is_snp14(self):
        self.assertEqual(self.case1.is_snp(min_count=3, min_freq=0.05,
                                           A=1,
                                           C=99,
                                           G=0,
                                           T=0,
                                           cov=100), (True, False) )

if __name__ == '__main__':
    unittest.main()
    
    
    