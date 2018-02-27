# -*- coding: utf-8 -*-

import unittest
import singleSL


class TestReaction(unittest.TestCase):
    '''Test class for asserting synthetic lethal reactions'''
    def test_singleSL(self):
        self.assertEqual(len(singleSL.singleSL('../Models/iAF1260.xml',
                                               0.01,
                                               '../Models/iAF1260_elimination_list.xml',
                                               'ATPM',
                                               'glpk')), 265)
        self.assertEqual(len(singleSL.singleSL('../Models/iNJ661.xml',
                                               0.01,
                                               '../Models/iNJ661_elimination_list.xml',
                                               'ATPM',
                                               'glpk')), 309)
        self.assertEqual(len(singleSL.singleSL('../Models/STM_v1.0.xml',
                                               0.01,
                                               '../Models/STM_v1.0_elimination_list.xml',
                                               'ATPM',
                                               'glpk')), 332)

if __name__ == '__main__':
    unittest.main()
