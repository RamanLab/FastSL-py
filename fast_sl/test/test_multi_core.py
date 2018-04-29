# -*- coding: utf-8 -*-

from __future__ import absolute_import

from os.path import abspath, dirname, join

import cobra

from fast_sl.multi_core.parallel_rxns import (
                                              parallel_single_sl,
                                              parallel_double_sl)
from fast_sl.multi_core.parallel_genes import (
                                               parallel_single_sl_genes,
                                               parallel_double_sl_genes)


model_dir_path = abspath(join(dirname(abspath(__file__)), "models", "stock",
                              "e_coli_core", "e_coli_core.xml"))
model = cobra.io.read_sbml_model(model_dir_path)
elilist = model.exchanges


def test_parallel_single_rxns():
    assert len(parallel_single_sl(model, 0.01, elilist, 'glpk_exact')) == 14


def test_parallel_double_rxns():
    assert len(parallel_double_sl(model, 0.01, elilist, 'glpk_exact')[1]) == 88


def test_parallel_single_genes():
    assert len(parallel_single_sl_genes(model, 0.01, 'glpk_exact')) == 7


def test_parallel_double_genes():
    assert len(parallel_double_sl_genes(model, 0.01, 'glpk_exact')[1]) == 53
