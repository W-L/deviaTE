#!/usr/bin/env python3

import unittest
import pysam
from deviaTE import deviaTE_multiHSP as multiHSP



class Test_build_cigar(unittest.TestCase):

    def setUp(self):
        self.case1 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=399, ref_end=499,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 120), (0, 100), (4, 280)],
                                                al_start=120, al_end=220,
                                                ref_start=599, ref_end=699,
                                                orig_container='sam_format')])
    
    
        self.case2 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=399, ref_end=499,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 120), (0, 100), (4, 280)],
                                                al_start=120, al_end=220,
                                                ref_start=499, ref_end=599,
                                                orig_container='sam_format')])


        self.case3 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=399, ref_end=499,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 120), (0, 100), (4, 280)],
                                                al_start=120, al_end=220,
                                                ref_start=494, ref_end=594,
                                                orig_container='sam_format')])
        
        
        self.case4 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=349, ref_end=449,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 110), (0, 100), (4, 290)],
                                                al_start=110, al_end=210,
                                                ref_start=549, ref_end=649,
                                                orig_container='sam_format')])
        
        
        self.case5 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=399, ref_end=499,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 110), (0, 100), (4, 290)],
                                                al_start=110, al_end=210,
                                                ref_start=499, ref_end=599,
                                                orig_container='sam_format')])
        
        
        self.case6 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=399, ref_end=499,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 110), (0, 100), (4, 290)],
                                                al_start=110, al_end=210,
                                                ref_start=494, ref_end=594,
                                                orig_container='sam_format')])
        
        
        self.case7 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=349, ref_end=449,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 105), (0, 100), (4, 295)],
                                                al_start=105, al_end=205,
                                                ref_start=549, ref_end=649,
                                                orig_container='sam_format')])
        
        
        self.case8 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=399, ref_end=499,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 105), (0, 100), (4, 295)],
                                                al_start=105, al_end=205,
                                                ref_start=499, ref_end=599,
                                                orig_container='sam_format')])
    
    
        self.case9 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(4, 10), (0, 100), (4, 390)],
                                                al_start=10, al_end=110,
                                                ref_start=399, ref_end=499,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 105), (0, 100), (4, 295)],
                                                al_start=105, al_end=205,
                                                ref_start=494, ref_end=594,
                                                orig_container='sam_format')])
        
        
        self.case10 = multiHSP.MAC(read_id='read_name', fam='unit_te',
                         hsp_list=[multiHSP.HSP(cigartuples=[(0, 100), (4, 400)],
                                                al_start=0, al_end=100,
                                                ref_start=999, ref_end=1099,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 120), (0, 100), (4, 280)],
                                                al_start=120, al_end=220,
                                                ref_start=1119, ref_end=1219,
                                                orig_container='sam_format'),
                                   multiHSP.HSP(cigartuples=[(4, 240), (0, 100), (4, 160)],
                                                al_start=240, al_end=340,
                                                ref_start=1239, ref_end=1339,
                                                orig_container='sam_format')])
        
        
    
    def test_build_cigar1(self):
        self.assertEqual(self.case1.build_cigar(), '10S100M10I100N100M280S')
        
    def test_build_cigar2(self):
        self.assertEqual(self.case2.build_cigar(), '10S100M10I100M280S')

    def test_build_cigar3(self):
        self.assertEqual(self.case3.build_cigar(), '10S100M10I95M285S')
        
    def test_build_cigar4(self):
        self.assertEqual(self.case4.build_cigar(), '10S100M100N100M290S')

    def test_build_cigar5(self):
        self.assertEqual(self.case5.build_cigar(), '10S200M290S')

    def test_build_cigar6(self):
        self.assertEqual(self.case6.build_cigar(), '10S195M295S')
        
    def test_build_cigar7(self):
        self.assertEqual(self.case7.build_cigar(), '10S100M105N95M295S')
        
    def test_build_cigar8(self):
        self.assertEqual(self.case8.build_cigar(), '10S100M5N95M295S')

    def test_build_cigar9(self):
        self.assertEqual(self.case9.build_cigar(), '10S195M295S')
        
    def test_build_cigar10(self):
        self.assertEqual(self.case10.build_cigar(), '100M20I20N100M20I20N100M160S')
        

if __name__ == '__main__':
    unittest.main()
