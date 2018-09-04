# -*- coding: utf-8 -*-

'''Contains pytest fixtures.'''

import pytest
from tests import ELILIST, MODEL


@pytest.fixture(scope="session")
def model():
    return MODEL


@pytest.fixture(scope="session")
def elilist():
    return ELILIST
