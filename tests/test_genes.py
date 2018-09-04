# -*- coding: utf-8 -*-

'''Contains test functions for fastsl genes.'''

import pytest

from fastsl.genes import single_genes, double_genes
from fastsl.parallel_genes import (
    parallel_single_genes,
    parallel_double_genes)


@pytest.mark.parametrize("mode_function", [single_genes,
                                           parallel_single_genes])
def test_single_genes(model, mode_function):
    '''Test single synthetic lethal genes.'''
    assert len(mode_function(model, 0.01, solver='glpk_exact')) == 7


@pytest.mark.parametrize("mode_function", [double_genes,
                                           parallel_double_genes])
def test_double_genes(model, mode_function):
    '''Test double synthetic lethal genes.'''
    assert len(mode_function(model, 0.01, solver='glpk_exact')[1]) == 53
