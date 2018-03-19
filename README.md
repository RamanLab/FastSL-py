# FastSL-py

This is the Python implementation of FastSL. CPLEX and Gurobi are not included.

## Basic requirement(s): 

    - Python 3.6.1 (conda distribution recommended)
    - Python 3.5.x for IBM CPLEX (conda distribution recommended)

## Installation:

- Check if your system has conda package manager installed before proceeding. If you don't have it, please install it from <https://www.continuum.io/downloads> by selecting desired Python distribution.

### Linux / macOS:
 
- In the terminal, run `sh install.sh`

### Windows:

- Open a terminal and run `conda env create -f environment.yml -n FastSL-py`.
 
- In the terminal, run `activate FastSL-py`.

## Usage:

- **Running FastSL-py on E.coli iAF1260**:
    In a terminal, run `python fast_sl.py 'Models/iAF1260.xml' --elilist 'Models/iAF1260_elimination_list.xml'`. 

- **Running FastSL-py genes on E.coli iAF1260**:
    Use the `--genes 1` flag.

- **Running the parallel version**:
    Use the `--parallel 1` flag.

- **Running other models not provided in the Models directory**:
    Use the `--gen_elilist True` flag.

- **For help**:
    Use the `-h` flag.
