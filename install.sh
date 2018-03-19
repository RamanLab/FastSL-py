#!/usr/bin/sh

# Initiate output directory:
if [ ! -d "Results" ]; then
  mkdir Results
fi

echo "Results directory created."

# Download dependencies:
conda env create -f environment.yml -n FastSL-py
echo "conda environment created."

# Initiating conda env:
echo "Activating conda environment."
source activate FastSL-py

