# FastSL-py

This is the Python implementation of FastSL.

## Basic requirement(s): 

    - Python 3.5 for IBM CPLEX
    - Python 3.6.1 (conda distribution recommended)

## Setting up:

1. **Installing conda**:
    Check if your system has conda package manager installed before proceeding. If you don't have it, please install it from <https://www.continuum.io/downloads> by selecting desired Python distribution.

2. **Creating conda environment**:
    Open a terminal and run `conda env create -f environment.yml -n FastSL-py`.
    
3. **Activating conda environment**:
    In the terminal, run `source activate FastSL-py`.

## Running the program:

1. **Running FastSL-py on E.coli iAF1260**:
    To run the program, navigate to `Rxns` directory and run `python fastSL.py '../Models/iAF1260.xml' --eliList '../Models/iAF1260_elimination_list.xml'`. 
    
2. **Running other models not provided in the Models directory**:
    To generate the elimination list, run `gen_eliList_xml.py`. After the elimination list generation, follow the same procedure as one given above.

