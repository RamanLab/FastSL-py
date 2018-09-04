# -*- coding: utf-8 -*-


'''Contains test functions for fastsl.'''

from fastsl.reactions import single_reactions, double_reactions
from fastsl.genes import single_genes, double_genes


def test_single_reactions(model, elilist):
    '''Test single synthetic lethal reactions.'''
    assert len(single_reactions(model, elilist, solver='glpk_exact')) == 14


def test_double_reactions(model, elilist):
    '''Test double synthetic lethal reactions.'''
    assert len(double_reactions(model, elilist, solver='glpk_exact')[1]) == 88


def test_single_genes(model):
    '''Test single lethal genes.'''
    assert len(single_genes(model, 0.01, solver='glpk_exact')) == 7


def test_double_genes(model):
    '''Test double lethal genes.'''
    assert len(double_genes(model, 0.01, solver='glpk_exact')[1]) == 53
