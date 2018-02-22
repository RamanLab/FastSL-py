# -*- coding: utf-8 -*-

import numpy as np
from tqdm import tqdm
import math
from singleSL import singleSL


def tripleSL(model, cutoff, eliList, atpm, solver):
    '''
    Analysis for triple lethal reactions
    '''

    # Wildtype FBA solution
    model.solver = solver
    solWT = model.optimize()
    grWT = solWT.objective_value
    Jnz = np.flatnonzero(solWT.fluxes)
    # If a list of reactions for which are eliminated for lethality is given
    # often exchange reactions are not considered

    # Index of reactions not considered for lethality analysis
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList]

    Jnz_copy = np.setdiff1d(Jnz, eliIdx)  # Jnz

    Jsl = singleSL(model=model, cutoff=cutoff, eliList=eliList, atpm=atpm,
                   solver=solver)
    Jsl_idx = [model.reactions.index(Jsl_id) for Jsl_id in Jsl]

    # Eliminate Single lethal reaction deletions for enumeration of higher
    # order lethals
    Jnz_copy = np.setdiff1d(Jnz_copy, Jsl_idx)

    Jdl_idx = []
    Jtl_idx = []

    for delIdx_i in tqdm(Jnz_copy, desc="Identifying Jdl and Jtl -\
            Part 1 of 2"):
        with model:
            model.reactions[delIdx_i].knock_out()
            # It can't be a single lethal so we can proceed further
            solKO_i = model.optimize()
            Jnz_i_before_primary_filtering = np.flatnonzero(solKO_i.fluxes)
            Jnz_i_before_secondary_filtering =\
                                                np.setdiff1d(Jnz_i_before_primary_filtering, Jnz)

            # Eliminate Exchange and ATP Maintenance reactions
            Jnz_i = np.setdiff1d(Jnz_i_before_secondary_filtering, eliIdx)

            for delIdx_j in Jnz_i:
                with model:
                    model.reactions[delIdx_j].knock_out()
                    solKO_ij = model.optimize()
                    if solKO_ij.objective_value < cutoff * grWT and \
                            solKO_ij.status != 'infeasible':
                        Jdl_idx.append([int(delIdx_i), int(delIdx_j)])
                    elif solKO_ij.status == 'infeasible':
                        solKO_ij = model.optimize()
                        if solKO_ij.objective_value < cutoff * grWT or \
                                math.isnan(solKO_ij.objective_value) is True:
                            Jdl_idx.append([int(delIdx_i), int(delIdx_j)])
                            # continue

                        Jnz_ij_before_primary_filtering = np.flatnonzero(solKO_ij.fluxes)
                        Jnz_ij_before_secondary_filtering = np.setdiff1d(Jnz_ij_before_primary_filtering,Jnz)

                        Jnz_ij = np.setdiff1d(Jnz_ij_before_secondary_filtering,eliIdx) # Eliminate Exchange and ATPM reactions

                        for delIdx_k in Jnz_ij:
                            with model:
                                model.reactions[delIdx_k].knock_out()
                                solKO_ijk = model.optimize()

                                if solKO_ijk.objective_value < cutoff * grWT or math.isnan(solKO_ijk.objective_value) == True:
                                    Jtl_idx.append([int(delIdx_i),int(delIdx_j),int(delIdx_k)])

    for delIdx_i in tqdm(Jnz_copy, desc='Identifying Jdl and Jtl -\
            Part 2 of 2'):
        for delIdx_j in Jnz_copy:
            if np.where(Jnz_copy == delIdx_j) < \
                    np.where(Jnz_copy == delIdx_i):
                with model:
                    model.reactions[delIdx_i].knock_out()
                    model.reactions[delIdx_j].knock_out()
                    solKO_ij = model.optimize()
                    if solKO_ij.objective_value < cutoff * grWT and solKO_ij.status != 'infeasible':
                        Jdl_idx.append([int(delIdx_i), int(delIdx_j)])
                    elif solKO_ij.status == 'infeasible':
                        solKO_ij = model.optimize()
                        if solKO_ij.objective_value < cutoff * grWT or math.isnan(solKO_ij.objective_value) == True:
                            Jdl_idx.append([int(delIdx_i), int(delIdx_j)])
                            # continue

                        Jnz_ij_before_primary_filtering =\
                                                           np.flatnonzero(solKO_ij.fluxes)
                        Jnz_ij_before_secondary_filtering =\
                                                            np.setdiff1d(Jnz_ij_before_primary_filtering,Jnz)
                        Jnz_ij = \
                                np.setdiff1d(Jnz_ij_before_secondary_filtering,
                                        eliIdx)  # Eliminate Exchange and ATPM reactions

                        for delIdx_k in Jnz_ij:
                            with model:
                                solKO_ijk = model.optimize()
                                if solKO_ijk.objective_value < cutoff * grWT or math.isnan(solKO_ijk.objective_value) == True:
                                    Jtl_idx.append([int(delIdx_i),int(delIdx_j),int(delIdx_k)])

                        for delIdx_k in Jnz_copy:
                            with model:
                                if np.where(Jnz_copy == delIdx_k) < \
                                        np.where(Jnz_copy == delIdx_j):
                                    solKO_ijk = model.optimize()
                                    if solKO_ijk.objective_value < \
                                            cutoff * grWT or \
                                            math.isnan(solKO_ijk.objective_value) == True:
                                        Jtl_idx.append([int(delIdx_i),int(delIdx_j),int(delIdx_k)])

                                else:
                                    break
            else:
                break

    # Eliminate double lethal reaction deletions in triple lethal reacutions
    # Jdl_idx = np.array(Jdl_idx)
    # Jtl_idx = np.array(Jtl_idx)

    # temporary = []
    # g = np.zeros(Jdl_idx.shape[0])
    # for delIdx_i in Jtl_idx:
    #     for delIdx_j in Jdl_idx:
    #         g[np.where(Jdl_idx == delIdx_j)[0][0]] = \
    #                 np.sum(np.in1d(delIdx_i, delIdx_j))
    #         if g[np.where(Jdl_idx == delIdx_j)[0][0]] >= 2:
    #             break
    #     if np.max(g) < 2:
    #         temporary.append(delIdx_i)

    # Jtl_idx_new = np.array(temporary)
    # Jtl_idx_final = np.unique(np.sort(Jtl_idx_new), axis=0)

    Jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx in Jdl_idx]
    Jtl = [model.reactions.get_by_any(rxn_triplet_idx) for rxn_triplet_idx in Jtl_idx]

    return (Jsl, Jdl, Jtl)
