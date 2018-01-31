#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import unittest
import singleSL

class TestReaction(unittest.TestCase):
    '''Test class for asserting synthetic lethal reactions'''
    def test_singleSL(self):
        self.assertEqual(len(singleSL.singleSL('../Models/iAF1260.xml',0.01,'../Models/iAF1260_elimination_list.xml','ATPM','glpk')),278)

if __name__ == '__main__':
    unittest.main()
