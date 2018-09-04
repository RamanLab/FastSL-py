# -*- coding: utf-8 -*-

'''Contains utilities for generating pytest fixtures.'''

from os.path import abspath, dirname, join

import cobra

MODEL_DIR_PATH = abspath(join(dirname(abspath(__file__)), "models", "stock",
                              "e_coli_core", "e_coli_core.xml"))
MODEL = cobra.io.read_sbml_model(MODEL_DIR_PATH)
ELILIST = MODEL.exchanges
