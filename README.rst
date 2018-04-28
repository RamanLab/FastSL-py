FastSL-py
=========

This is the Python implementation of FastSL, an efficient algorithm to
identify synthetic lethal gene/reaction sets in genome-scale metabolic
models.

This package is based on
`cobrapy <https://github.com/opencobra/cobrapy>`__ and provides a simple
command-line tool.

If you use FastSL-py in your work, please do cite:
doi:10.1093/bioinformatics/btv352

Basic requirement(s):
---------------------

::

    - Python 3.6 for Gurobi 8 (conda distribution recommended)
    - Python 3.5 for IBM CPLEX and Gurobi 7 (conda distribution recommended)

Installation:
-------------

Use pip to install from PyPI (recommended inside a virtual environment):

::

    pip install fast_sl

In case you downloaded the source code from GitHub:

::

    pip install -e .

Usage:
------

-  **Running FastSL-py on a COBRA model (.xml file)**: In a terminal,
   run
   ``python fast_sl.py '<model.xml>' --elilist '<model_elimination_list.xml>'``.

-  **Running FastSL-py genes on E.coli iAF1260**: Use the ``--genes 1``
   flag.

-  **Running the parallel version**: Use the ``--parallel 1`` flag.

-  **Generating elimination list**: Use the ``--gen_elilist True`` flag.

-  **For help**: Use the ``-h`` flag.

Note:
-----

CPLEX and Gurobi are not included. Both are available for free (for
academic purposes). All solvers are supported whose interfaces are
provided by `optlang <https://github.com/biosustain/optlang>`__.
