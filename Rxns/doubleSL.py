#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import numpy
from tqdm import tqdm
import math
from singleSL import singleSL

def doubleSL(model,cutoff,eliList,atpm,solver):
    '''
    Analysis for double lethal reactions
    '''
    model.solver = solver
    solWT = model.optimize()
    grWT = solWT.objective_value
    J = solWT.fluxes
    Jnz = numpy.flatnonzero(J)
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList]
    Jnz = numpy.setdiff1d(Jnz,eliIdx) # Jnz

    Jsl = singleSL(model=model,cutoff=cutoff,eliList=eliList,atpm=atpm,solver=solver)
    Jsl_idx = [model.reactions.index(Jsl_id) for Jsl_id in Jsl]

    Jnz_copy = numpy.setdiff1d(Jnz,Jsl_idx) # Jnz-Jsl

    Jdl_idx = []

    for delIdx_i in tqdm(Jnz_copy,desc="Identifying Jdl - Part 1 of 2"):
        with model:
            model.reactions[delIdx_i].knock_out()
            solKO_i = model.optimize()
            newnnz = numpy.flatnonzero(solKO_i.fluxes)
            Jnz_i = numpy.setdiff1d(newnnz,Jnz)
            Jnz_i = numpy.setdiff1d(Jnz_i,eliIdx)

            for delIdx_j in Jnz_i:
                with model:
                    model.reactions[delIdx_j].knock_out()
                    solKO_ij = model.slim_optimize()
                    if solKO_ij < cutoff * grWT or math.isnan(solKO_ij) == True:
                        Jdl_idx.append([int(delIdx_i),int(delIdx_j)])

    for delIdx_i in tqdm(Jnz_copy,desc="Identifying Jdl - Part 2 of 2"):
        for delIdx_j in Jnz_copy:
            if numpy.where(Jnz_copy==delIdx_j) < numpy.where(Jnz_copy==delIdx_i):
                with model:
                    model.reactions[delIdx_i].knock_out()
                    model.reactions[delIdx_j].knock_out()
                    solKO_ij = model.slim_optimize()
                    if solKO_ij < cutoff * grWT or math.isnan(solKO_ij) == True:
                        Jdl_idx.append([int(delIdx_i),int(delIdx_j)])

    #Jsl = model.reactions.get_by_any(Jsl)
    Jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx in Jdl_idx]
    return (Jsl,Jdl)
