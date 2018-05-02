FastSL-py
=========

.. image:: https://travis-ci.org/RamanLab/FastSL-py.svg?branch=master
    :target: https://travis-ci.org/RamanLab/FastSL-py
    :alt: Travis CI Status

.. image::  https://ci.appveyor.com/api/projects/status/d97plfb60bcic5ar?svg=true
    :target: https://ci.appveyor.com/project/synchon/fastsl-py
    :alt: AppVeyor CI Status

.. image:: https://readthedocs.org/projects/fastsl-py/badge/?version=latest
    :target: http://fastsl-py.readthedocs.io/?badge=latest
    :alt: Documentation Status


This is the Python implementation of FastSL, an efficient algorithm to
identify synthetic lethal gene/reaction sets in genome-scale metabolic
models.

This package is based on
`cobrapy <https://github.com/opencobra/cobrapy>`__ and provides a simple
command-line tool.

For documentation, please visit: `http://fastsl-py.readthedocs.io`

Basic requirement(s):
---------------------

::

    - Python 3.6 for Gurobi 8
    - Python 3.5 for IBM CPLEX and Gurobi 7

Installation:
-------------

Use pip to install from PyPI (recommended inside a virtual environment):

::

    pip install fastsl

Contribute:
-----------

- Issue Tracker: `<https://github.com/RamanLab/FastSL-py/issues>`

Support:
--------

If you are having issues, please let us know.
Contact us at: <fast-sl@ramanlab.groups.io>

License:
--------

The project is licensed under GPL v3 license.

Note:
-----

CPLEX and Gurobi are not included. Both are available for free (for
academic purposes). All solvers are supported whose interfaces are
provided by `optlang <https://github.com/biosustain/optlang>`__.
