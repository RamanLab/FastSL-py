#!/bin/zsh

# Takes a root directory which has the folders comprising of model and elimination list files
# and generates the lethals sets

ROOT_DIR=$1;

for SUB_DIR in "$ROOT_DIR"/*
do
    for MODEL_FILES in "$SUB_DIR"/*
    do
        if [[ "$MODEL_FILES" == *_elimination_list.xml ]]
        then
            elilist="$MODEL_FILES";
        else
            model="$MODEL_FILES"
        fi
    done
    
    # Change the parameter(s) as required
    
    # run reactions
    python fast_sl.py "$model" --elilist "$elilist" --solver $2 --parallel 1 --genes 0
    
    # run genes
    python fast_sl.py "$model" --solver $2 --parallel 1 --genes 1

done
