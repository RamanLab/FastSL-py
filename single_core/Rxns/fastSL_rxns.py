# -*- coding: utf-8 -*-

import math
from itertools import product

import numpy as np
from tqdm import tqdm


def singleSL(model, cutoff, eliList, solver):
    '''
    Analysis for single lethal reactions
    '''

    model.solver = solver  # Basic solver configuration
    solWT = model.optimize()  # Identify minNorm flux distribution
    grWT = solWT.objective_value

    # Indices of non-zero flux reactions including exchange/diffusion reactions
    Jnz_before_filtering = np.flatnonzero(solWT.fluxes)

    # Indices of exchange/diffusion reactions (not considered)
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList]

    Jnz = np.setdiff1d(Jnz_before_filtering, eliIdx)  # Jnz

    # Identify Single Lethal Reaction Deletions

    Jsl_idx = []

    for delIdx_i in tqdm(Jnz, desc="Identifying Jsl reactions"):
        with model:
            model.reactions[delIdx_i].knock_out()
            solKO_i = model.slim_optimize()
            if solKO_i < cutoff * grWT or math.isnan(solKO_i) is True:
                Jsl_idx.append(int(delIdx_i))

    # Indices -> Reactions
    Jsl = model.reactions.get_by_any(Jsl_idx)

    return Jsl


def doubleSL(model, cutoff, eliList, solver):
    '''
    Analysis for double lethal reactions
    '''

    model.solver = solver  # Basic solver configuration
    solWT = model.optimize()  # Identify minNorm flux distribution
    grWT = solWT.objective_value

    # Indices of non-zero flux reactions including exchange/diffusion reactions
    Jnz_before_filtering = np.flatnonzero(solWT.fluxes)

    # Indices of exchange/diffusion reactions (not considered)
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList]

    Jnz = np.setdiff1d(Jnz_before_filtering, eliIdx)  # Jnz

    # Identify single lethal reactions
    Jsl = singleSL(model,
                   cutoff,
                   eliList,
                   solver)

    # Indices of single lethal reacions
    Jsl_idx = [model.reactions.index(Jsl_id) for Jsl_id in Jsl]

    Jnz_copy = np.setdiff1d(Jnz, Jsl_idx)  # Jnz-Jsl

    # Makes rxn pairs to remove nested for-loops in phase 2
    Jnz_copy_phase_2 = [rxn_pair for rxn_pair in product(Jnz_copy, repeat=2)]

    # Identify double lethal reactions

    Jdl_idx = []

    # Phase 1:
    for delIdx_i in tqdm(Jnz_copy, desc="Identifying Jdl reactions: 1 of 2"):
        with model:
            model.reactions[delIdx_i].knock_out()
            solKO_i = model.optimize()
            newnnz = np.flatnonzero(solKO_i.fluxes)
            Jnz_i_before_filtering = np.setdiff1d(newnnz, Jnz)
            Jnz_i = np.setdiff1d(Jnz_i_before_filtering, eliIdx)

            for delIdx_j in Jnz_i:
                with model:
                    model.reactions[delIdx_j].knock_out()
                    solKO_ij = model.slim_optimize()
                    if solKO_ij < cutoff * grWT or \
                            math.isnan(solKO_ij) is True:
                        Jdl_idx.append([int(delIdx_i), int(delIdx_j)])

    # Phase 2:
    for delIdx_pair in tqdm(Jnz_copy_phase_2,
                            desc="Identifying Jdl reactions: 2 of 2"):
        delIdx_i, delIdx_j = delIdx_pair
        if np.where(Jnz_copy == delIdx_j) < np.where(Jnz_copy == delIdx_i):
            with model:
                model.reactions[delIdx_i].knock_out()
                model.reactions[delIdx_j].knock_out()
                solKO_ij = model.slim_optimize()
                if solKO_ij < cutoff * grWT or \
                        math.isnan(solKO_ij) is True:
                    Jdl_idx.append([int(delIdx_i), int(delIdx_j)])

    # Indices -> Reactions
    Jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in Jdl_idx]

    return (Jsl, Jdl)


def tripleSL(model, cutoff, eliList, solver):
    '''
    Analysis for triple lethal reactions
    '''

    model.solver = solver  # Basic solver configuration
    solWT = model.optimize()  # Identify minNorm flux distribution
    grWT = solWT.objective_value

    # Indices of non-zero flux reactions including exchange/diffusion reactions
    Jnz_before_filtering = np.flatnonzero(solWT.fluxes)

    # Indices of exchange/diffusion reactions (not considered)
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList]

    Jnz = np.setdiff1d(Jnz_before_filtering, eliIdx)  # Jnz

    # Identify single lethal reactions
    Jsl = singleSL(model,
                   cutoff,
                   eliList,
                   solver)

    # Indices of single lethal reacions
    Jsl_idx = [model.reactions.index(Jsl_id) for Jsl_id in Jsl]

    Jnz_copy = np.setdiff1d(Jnz, Jsl_idx)  # Jnz-Jsl

    # Makes rxn pairs to remove nested for-loops in phase 2
    Jnz_copy_phase_2 = [rxn_pair for rxn_pair in product(Jnz_copy, repeat=2)]

    # Identify double lethal reactions

    Jdl_idx = []
    Jtl_idx = []

    # Phase 1:
    for delIdx_i in tqdm(Jnz_copy,
                         desc="Identifying Jdl and Jtl reactions: 1 of 2"):
        with model:
            model.reactions[delIdx_i].knock_out()
            solKO_i = model.optimize()
            newnnz_i = np.flatnonzero(solKO_i.fluxes)
            Jnz_i_before_filtering = np.setdiff1d(newnnz_i, Jnz)
            Jnz_i = np.setdiff1d(Jnz_i_before_filtering, eliIdx)

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

                        newnnz_ij = np.flatnonzero(solKO_ij.fluxes)
                        Jnz_ij_before_filtering = np.setdiff1d(newnnz_ij, Jnz)
                        Jnz_ij = np.setdiff1d(Jnz_ij_before_filtering, eliIdx)

                        for delIdx_k in Jnz_ij:
                            with model:
                                model.reactions[delIdx_k].knock_out()
                                solKO_ijk = model.slim_optimize()
                                if solKO_ijk < cutoff * grWT or \
                                        math.isnan(solKO_ijk) is True:
                                    Jtl_idx.append([int(delIdx_i),
                                                    int(delIdx_j),
                                                    int(delIdx_k)])

    # Phase 2:
    for delIdx_pair in tqdm(Jnz_copy_phase_2,
                            desc="Identifying Jdl and Jtl reactions: 2 of 2"):
        delIdx_i, delIdx_j = delIdx_pair
        if np.where(Jnz_copy == delIdx_j) < np.where(Jnz_copy == delIdx_i):
            with model:
                model.reactions[delIdx_i].knock_out()
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

                    newnnz_ij = np.flatnonzero(solKO_ij.fluxes)
                    Jnz_ij_before_filtering = np.setdiff1d(newnnz_ij, Jnz)
                    Jnz_ij = np.setdiff1d(Jnz_ij_before_filtering, eliIdx)

                    for delIdx_k in Jnz_ij:
                        with model:
                            solKO_ijk = model.slim_optimize()
                            if solKO_ijk < cutoff * grWT or \
                                    math.isnan(solKO_ijk) is True:
                                Jtl_idx.append([int(delIdx_i),
                                                int(delIdx_j),
                                                int(delIdx_k)])

                    for delIdx_k in Jnz_copy:
                        with model:
                            if np.where(Jnz_copy == delIdx_k) < \
                                    np.where(Jnz_copy == delIdx_j):
                                solKO_ijk = model.slim_optimize()
                                if solKO_ijk < cutoff * grWT or \
                                        math.isnan(solKO_ijk) is True:
                                    Jtl_idx.append([int(delIdx_i),
                                                    int(delIdx_j),
                                                    int(delIdx_k)])

    # Eliminate double lethal reaction deletions in triple lethal reacutions
    Jdl_idx = np.array(Jdl_idx)
    Jtl_idx = np.array(Jtl_idx)

    temporary = []
    g = np.zeros(Jdl_idx.shape[0])
    for delIdx_i in Jtl_idx:
        for delIdx_j in Jdl_idx:
            g[np.where(Jdl_idx == delIdx_j)[0][0]] = \
                    np.sum(np.isin(delIdx_i, delIdx_j))
            # if g[np.where(Jdl_idx == delIdx_j)[0][0]] >= 2:
            #     continue
        if np.max(g) < 2:
            temporary.append(delIdx_i)

    Jtl_idx_new = np.array(temporary)
    Jtl_idx_final = np.unique(np.sort(Jtl_idx_new), axis=0).tolist()

    # Indices -> Reactions
    Jdl_idx = Jdl_idx.tolist()
    Jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in Jdl_idx]

    Jtl = [model.reactions.get_by_any(rxn_triplet_idx) for rxn_triplet_idx
           in Jtl_idx_final]

    return (Jsl, Jdl, Jtl_idx_final)
