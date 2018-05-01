# -*- coding: utf-8 -*-

from os.path import abspath, dirname, join

import cobra

from fastsl.rxns import single_sl, double_sl
from fastsl.genes import single_sl_genes, double_sl_genes

MODEL_DIR_PATH = abspath(join(dirname(abspath(__file__)), "models", "stock",
                              "e_coli_core", "e_coli_core.xml"))
MODEL = cobra.io.read_sbml_model(MODEL_DIR_PATH)
ELILIST = MODEL.exchanges


def test_single_rxns():
    '''Test function for single reaction deletions.'''
    assert len(single_sl(MODEL, 0.01, ELILIST, 'glpk_exact')) == 14


def test_double_rxns():
    '''Test function for double reaction deletions.'''
    assert len(double_sl(MODEL, 0.01, ELILIST, 'glpk_exact')[1]) == 88


def test_single_genes():
    '''Test function for single gene deletions.'''
    assert len(single_sl_genes(MODEL, 0.01, 'glpk_exact')) == 7


def test_double_genes():
    '''Test function for double gene deletions.'''
    assert len(double_sl_genes(MODEL, 0.01, 'glpk_exact')[1]) == 53
