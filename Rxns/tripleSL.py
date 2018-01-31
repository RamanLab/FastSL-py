#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import numpy
from tqdm import tqdm
import math
from singleSL import singleSL

def tripleSL(model,cutoff,eliList,atpm,solver):
    '''
    Analysis for triple lethal reactions
    '''

    # Wildtype FBA solution
    model.solver = solver
    solWT = model.optimize()
    grWT = solWT.objective_value
    J = solWT.fluxes
    Jnz = numpy.flatnonzero(J)
    # If a list of reactions for which are eliminated for lethality is given often exchange reactions are not considered
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList] # Index of reactions not considered for lethality analysis
    Jnz_copy = numpy.setdiff1d(Jnz,eliIdx) # Jnz

    Jsl = singleSL(model=model,cutoff=cutoff,eliList=eliList,atpm=atpm,solver=solver)
    Jsl_idx = [model.reactions.index(Jsl_id) for Jsl_id in Jsl]
    Jnz_copy = numpy.setdiff1d(Jnz_copy,Jsl_idx) # Eliminate Single lethal reaction deletions for enumeration of higher order lethals

    Jdl_idx = []
    Jtl_idx = []

    for delIdx_i in tqdm(Jnz_copy,desc="Identifying Jdl and Jtl - Part 1 of 2"):
        with model:
            model.reactions[delIdx_i].knock_out()
            solKO_i = model.optimize() # It can't be a single lethal so we can proceed further
            Jnz_i = numpy.flatnonzero(solKO_i.fluxes)
            Jnz_i = numpy.setdiff1d(Jnz_i,Jnz)
            Jnz_i = numpy.setdiff1d(Jnz,eliIdx) # Eliminate Exchange and ATP Maintenance reactions

            for delIdx_j in Jnz_i:
                with model:
                    model.reactions[delIdx_j].knock_out()
                    solKO_ij = model.optimize()
                    if solKO_ij.objective_value < cutoff * grWT and solKO_ij.status != 'infeasible':
                        Jdl_idx.append([int(delIdx_i),int(delIdx_j)])
                    elif solKO_ij.status == 'infeasible':
                        solKO_ij = model.optimize()
                        if solKO_ij.objective_value < cutoff * grWT or math.isnan(solKO_ij.objective_value) == True:
                            Jdl_idx.append([int(delIdx_i),int(delIdx_j)])

                        Jnz_ij = numpy.flatnonzero(solKO_ij.fluxes)
                        Jnz_ij = numpy.setdiff1d(Jnz_ij,Jnz)

                        Jnz_ij = numpy.setdiff1d(Jnz_ij,eliIdx) # Eliminate Exchange and ATPM reactions

                        for delIdx_k in Jnz_ij:
                            with model:
                                model.reactions[delIdx_k].knock_out()
                                solKO_ijk = model.optimize()

                                if solKO_ijk.objective_value < cutoff * grWT or math.isnan(solKO_ijk.objective_value) == True:
                                    Jtl_idx.append([int(delIdx_i),int(delIdx_j),int(delIdx_k)])

    for delIdx_i in tqdm(Jnz_copy,desc='Identifying Jdl and Jtl = Part 2 of 2'):
        for delIdx_j in Jnz_copy:
            if numpy.where(Jnz_copy==delIdx_j) < numpy.where(Jnz_copy==delIdx_i):
                with model:
                    model.reactions[delIdx_i].knock_out()
                    model.reactions[delIdx_j].knock_out()
                    solKO_ij = model.optimize()
                    if solKO_ij.objective_value < cutoff * grWT and solKO_ij.status != 'infeasible':
                        Jdl_idx.append([int(delIdx_i),int(delIdx_j)])
                    elif solKO_ij.status == 'infeasible':
                        solKO_ij = model.optimize()
                        if solKO_ij.objective_value < cutoff * grWT or math.isnan(solKO_ij.objective_value) == True:
                            Jdl_idx.append([int(delIdx_i),int(delIdx_j)])

                        Jnz_ij = numpy.flatnonzero(solKO_ij.fluxes)
                        Jnz_ij = numpy.setdiff1d(Jnz_ij,Jnz)
                        Jnz_ij = numpy.setdiff1d(Jnz_ij,eliIdx) # Eliminate Exchange and ATPM reactions

                        for delIdx_k in Jnz_ij:
                            with model:
                                solKO_ijk = model.optimize()
                                if solKO_ijk.objective_value < cutoff * grWT or math.isnan(solKO_ijk.objective_value) == True:
                                    Jtl_idx.append([int(delIdx_i),int(delIdx_j),int(delIdx_k)])

                        for delIdx_k in Jnz_copy:
                            with model:
                                if numpy.where(Jnz_copy==delIdx_k) < numpy.where(Jnz_copy==delIdx_j):
                                    solKO_ijk = model.optimize()
                                    if solKO_ijk.objective_value < cutoff * grWT or math.isnan(solKO_ijk.objective_value) == True:
                                        Jtl_idx.append([int(delIdx_i),int(delIdx_j),int(delIdx_k)])

    # Eliminate double lethal reaction deletions in triple lethal reacutions
    Jdl_idx = numpy.array(Jdl_idx)
    Jtl_idx = numpy.array(Jtl_idx)

    temporary = []
    g = numpy.zeros((1,len(Jdl_idx)))
    for delIdx_i in Jtl_idx:
        for delIdx_j in Jdl_idx:
            g[numpy.where(delIdx_j)] = numpy.sum(numpy.in1d(Jtl_idx[numpy.where(delIdx_i),:],Jdl_idx[numpy.where(delIdx_j),:]))
            if g[numpy.where(delIdx_j)] >= 2:
                break
        if numpy.max(g) < 2:
            temporary.append(delIdx_i.tolist())
    
    Jtl_idx = numpy.array(temporary)
    Jtl_idx = numpy.unique(numpy.sort(Jtl_idx), axis=0)

    Jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx in Jdl_idx]
    Jtl = [model.reactions.get_by_any(rxn_triplet_idx) for rxn_triplet_idx in Jtl_idx]
    
    return (Jsl,Jdl,Jtl)
