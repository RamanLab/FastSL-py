# FastSL-py

This is the Python implementation of FastSL.

## Basic requirement(s): 

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
    To run the program, navigate to `single_core` directory and run `python fast_sl.py '../Models/iAF1260.xml' --elilist '../Models/iAF1260_elimination_list.xml'`. 

2. **Running FastSL-py genes on E.coli iAF1260**:
    To run the program, navigate to `single_core` directory and run  `python fast_sl.py '../Models/iAF1260.xml' --genes 1`.

3. **Running the parallel version**:
    Please navigate to the `multi_core` directory and follow the same command as described above but replace `fast_sl.py` with `parallel_fast_sl.py`.

4. **Running other models not provided in the Models directory**:
    To generate the elimination list, run `python gen_eliList_xml.py <path/to/model>`. After the elimination list generation, follow the same procedure as given above.

5. **For help**:
    Run the program with `-h` flag for detailed guide.
