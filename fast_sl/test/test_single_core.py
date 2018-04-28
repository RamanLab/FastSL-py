# -*- coding: utf-8 -*-

from __future__ import absolute_import


import cobra

from fast_sl.single_core.rxns import single_sl, double_sl
from fast_sl.single_core.genes import (
                                       single_sl_genes,
                                       double_sl_genes)


model = cobra.io.read_sbml_model('models/stock/e_coli_core/e_coli_core.xml')
elilist = model.exchanges


def test_single_rxns():
    assert len(single_sl(model, 0.01, elilist, 'glpk_exact')) == 14


def test_double_rxns():
    assert len(double_sl(model, 0.01, elilist, 'glpk_exact')[1]) == 88


def test_single_genes():
    assert len(single_sl_genes(model, 0.01, 'glpk_exact')) == 7


def test_double_genes():
    assert len(double_sl_genes(model, 0.01, 'glpk_exact')[1]) == 53
