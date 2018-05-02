Documentation for FastSL-py
===========================

.. toctree::
   :maxdepth: 1

FastSL-py is the Python implementation for the `FastSL <https://github.com/RamanLab/FastSL>`__ algorithm. This algorithm uses critical considerations which allows it to be highly efficient. This package provides a command-line tool which makes it easy to use the algorithm to obtain results for synthetic lethality analysis.

Installation
============

The recommended way is to install from PyPI (virtual environment recommended):

::

   >>> pip install fastsl

Usage
=====

Every operation is logged in **fast-sl.log** file which will be generated the first time the tool is used in the current directory.

Identifying synthetic lethal reactions
--------------------------------------

::

    >>> fast-sl <model-path> --elilist <eliminationlist-path>

This will generate a directory named **fast-sl-results** (if not generated already) and in that the **reactions** directory for the specific model will be generated containing the *.csv* files.
The elimination list is optional as if not provided it will currently consider the ATP maintenance reaction for elimination.

Identifying synthetic lethal genes
----------------------------------

::

    >>> fast-sl <model-path> --genes

This will generate a **genes** directory under the specific model directory. The **genes** directory will contain the required *.csv* files.

Using parallel mode
-------------------

For reactions:

::

    >>> fast-sl <model-path> --parallel

For genes:

::

    fast-sl <model-path> --parallel --genes

It performs better than the serial one and should be considered for simulating large models or bulk models.
This feature currently restricts the use to a max of 4 processes for now due to the bias in result found due to memory issues.
This feature is currently being worked on and will be capable of utilizing more processes in the future.

Generating elimination list
---------------------------

::

    >>> fast-sl <model-path> --gen-elilist

This will generate **<model-name>_elimination_list.xml** file in the current directory.

Specifying the cutoff value
-----------------------------------------------------------------------------------------------

::

    >>> fast-sl <model-path> --cutoff 0.1

This is used to set the fraction of growth rate to be considered for cutoff value of synthetic lethality. By default, 1% is considered as the cutoff value, but can be modified by passing in a different value to the option.

Specifying the order of synthetic lethals
-----------------------------------------

::

    >>> fast-sl <model-path> --order 1

By default, it has an order of 2 but can be changed to 1. In future, it will support orders 3 and 4 too.

Specifying the ID of ATP maintenance reaction
---------------------------------------------

::

    >>> fast-sl <model-path> --atpm 'ATPM'

By default, it takes ATPM as the ATP maintenance reaction ID for the model.

Specifying the LP solver to be used
-----------------------------------

::

    >>> fast-sl <model-path> --solver 'gurobi'

This allows to choose the desired LP solver which are supported. By default, it is 'glpk_exact'.
