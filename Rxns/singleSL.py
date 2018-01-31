#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import numpy as np
from tqdm import tqdm
import math
#import cobra.flux_analysis as cflux

def singleSL(model,cutoff,eliList,atpm,solver):
    '''
    Analysis for single lethal reactions
    '''

    # Step 1: Identify single lethal reactions
    # Identify minNorm flux distribution
    model.solver = solver
    solWT = model.optimize()
    grWT = solWT.objective_value
    #Jnz = solWT.fluxes[solWT.fluxes != 0].drop(eliList, errors='ignore').index
    #Jsl = cflux.single_reaction_deletion(model, Jnz)['flux'].where(lambda x : x < grWT * 0.01).dropna()
    Jnz = np.flatnonzero(solWT.fluxes)
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList] # Index of reactions not considered for lethality analysis
    Jnz = np.setdiff1d(Jnz,eliIdx) # Jnz

    # Identify Single Lethal Reaction Deletions...

    Jsl_idx = []

    for delIdx_i in tqdm(Jnz,desc="Identifying Jsl"):
        with model:
            model.reactions[delIdx_i].knock_out()
            solKO_i = model.slim_optimize()
            if solKO_i < cutoff * grWT or math.isnan(solKO_i) == True:
                Jsl_idx.append(int(delIdx_i))

    Jsl = model.reactions.get_by_any(Jsl_idx)
    return Jsl
