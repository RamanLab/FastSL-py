# -*- coding: utf-8 -*-

import math
from os import cpu_count
from itertools import product

import numpy as np
from joblib import Parallel, delayed

from singleSL import singleSL


def _double_lethal_reactions_phase_1_helper(model, cutoff, grWT, delIdx_i,
                                            delIdx_j):
    with model:
        model.reactions[delIdx_i].knock_out()
        model.reactions[delIdx_j].knock_out()
        solKO_ij = model.slim_optimize()
        if solKO_ij < cutoff * grWT or \
                math.isnan(solKO_ij) is True:
            return [int(delIdx_i), int(delIdx_j)]
        else:
            return None


def _double_lethal_reactions_phase_1(model, cutoff, grWT, Jnz, eliIdx,
                                     delIdx_i):
    with model:
        model.reactions[delIdx_i].knock_out()
        solKO_i = model.optimize()
        newnnz = np.flatnonzero(solKO_i.fluxes)
        Jnz_i_before_filtering = np.setdiff1d(newnnz, Jnz)
        Jnz_i = np.setdiff1d(Jnz_i_before_filtering, eliIdx)

    return list(
                filter(
                       lambda rxn_pair_idx: rxn_pair_idx is not None,
                       Parallel(n_jobs=1,
                                backend='multiprocessing',
                                verbose=0,
                                batch_size='auto')(
                                delayed(_double_lethal_reactions_phase_1_helper)
                                (model,
                                 cutoff,
                                 grWT,
                                 delIdx_i,
                                 delIdx_j) for delIdx_j in Jnz_i)))


def _double_lethal_reactions_phase_2(model, cutoff, grWT, Jnz_copy,
                                     delIdx_pair):
    delIdx_i, delIdx_j = delIdx_pair
    if np.where(Jnz_copy == delIdx_j) < np.where(Jnz_copy == delIdx_i):
        with model:
            model.reactions[delIdx_i].knock_out()
            model.reactions[delIdx_j].knock_out()
            solKO_ij = model.slim_optimize()
            if solKO_ij < cutoff * grWT or \
                    math.isnan(solKO_ij) is True:
                return [int(delIdx_i), int(delIdx_j)]
            else:
                return None


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

    Jnz_copy_phase_2 = [rxn_pair for rxn_pair in product(Jnz_copy, repeat=2)]

    # Identify Double Lethal Reaction Deletions

    # Phase 1
    chunk_size_phase_1 = Jnz_copy.shape[0] // cpu_count()
    Jdl_idx_1 = list(
                     filter(
                            lambda rxn_pair_idx: rxn_pair_idx is not None,
                            Parallel(n_jobs=cpu_count(),
                                     backend='multiprocessing',
                                     verbose=3,
                                     batch_size=chunk_size_phase_1)(
                                     delayed(_double_lethal_reactions_phase_1)
                                     (model,
                                      cutoff,
                                      grWT,
                                      Jnz,
                                      eliIdx,
                                      delIdx_i) for delIdx_i in Jnz_copy)))

    # Phase 2
    chunk_size_phase_2 = len(Jnz_copy_phase_2) // cpu_count()
    Jdl_idx_2 = list(
                     filter(
                            lambda rxn_idx: rxn_idx is not None,
                            Parallel(n_jobs=cpu_count(),
                                     backend='threading',
                                     verbose=3,
                                     batch_size=chunk_size_phase_2)(
                                     delayed(_double_lethal_reactions_phase_2)
                                     (model,
                                      cutoff,
                                      grWT,
                                      Jnz_copy,
                                      delIdx_pair) for delIdx_pair in Jnz_copy_phase_2)))

    Jdl = [inner_list for outer_list in list(filter(None, Jdl_idx_1))
           for inner_list in outer_list] + Jdl_idx_2

    # Indices -> Reactions
    Jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in Jdl_idx]

    return (Jsl, Jdl)
