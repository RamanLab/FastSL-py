# -*- coding: utf-8 -*-

'''Contains test functions for fastsl reactions.'''

import pytest

from fastsl.reactions import single_reactions, double_reactions
from fastsl.parallel_reactions import (
    parallel_single_reactions,
    parallel_double_reactions)


@pytest.mark.parametrize("mode_function", [single_reactions,
                                           parallel_single_reactions])
def test_single_reactions(model, elilist, mode_function):
    '''Test single synthetic lethal reactions.'''
    assert len(mode_function(model, elilist, solver='glpk_exact')) == 14


@pytest.mark.parametrize("mode_function", [double_reactions,
                                           parallel_double_reactions])
def test_double_reactions(model, elilist, mode_function):
    '''Test double synthetic lethal reactions.'''
    assert len(mode_function(model, elilist, solver='glpk_exact')[1]) == 88
