# -*- coding: utf-8 -*-

'''Contains functions to identify synthetic lethal reactions.'''

import math
from itertools import product
from os import cpu_count

import numpy as np
from joblib import Parallel, delayed


def _single_reactions_helper(model, cutoff, gr_wt, i):
    with model:
        model.reactions[i].knock_out()
        sol_ko_i = model.slim_optimize()
        if sol_ko_i < cutoff * gr_wt or math.isnan(sol_ko_i) is True:
            return int(i)


def _double_reactions_phase_1_helper(model, cutoff, gr_wt, i, j):
    with model:
        model.reactions[i].knock_out()
        model.reactions[j].knock_out()
        sol_ko_ij = model.slim_optimize()
        if sol_ko_ij < cutoff * gr_wt or \
                math.isnan(sol_ko_ij) is True:
            return [int(i), int(j)]


def _double_reactions_phase_1(model, cutoff, gr_wt, tolerance, jnz,
                              eli_idx, i):
    with model:
        model.reactions[i].knock_out()
        sol_ko_i = model.optimize()
        newnnz = np.where(sol_ko_i.fluxes.abs() > tolerance)[0]
        jnz_i_before_filtering = np.setdiff1d(newnnz, jnz)
        jnz_i = np.setdiff1d(jnz_i_before_filtering, eli_idx)

    return list(
                filter(
                       lambda rxn_pair_idx: rxn_pair_idx is not None,
                       Parallel(n_jobs=1,
                                backend='multiprocessing',
                                verbose=0,
                                batch_size='auto')(
                                delayed(
                                    _double_reactions_phase_1_helper)
                                (model, cutoff, gr_wt, i, j)
                                    for j in jnz_i)))


def _double_reactions_phase_2(model, cutoff, gr_wt, jnz_minus_jsl, pair):
    i, j = pair
    if np.where(jnz_minus_jsl == j) < np.where(jnz_minus_jsl == i):
        with model:
            model.reactions[i].knock_out()
            model.reactions[j].knock_out()
            sol_ko_ij = model.slim_optimize()
            if sol_ko_ij < cutoff * gr_wt or \
                    math.isnan(sol_ko_ij) is True:
                return [int(i), int(j)]


def parallel_single_reactions(model, elilist, cutoff=0.01, **kwargs):
    '''
    Analyze single lethal reactions.

    Parameters
    ----------
    model: cobra.Model
        The model to perform operation on.
    elilist: list
        The list to be used for discarding the reactions.
    cutoff: float, optional (default 0.01)
        The value to be used as growth cutoff for deciding
        synthetic lethality.
    **kwargs

    Returns
    -------
    list of single synthetic lethal reaction IDs

    '''
    solver = 'glpk'
    tolerance = model.solver.configuration.tolerances.optimality

    model.solver = solver
    sol_wt = model.optimize()
    gr_wt = sol_wt.objective_value


    # Indices of non-zero flux reactions including exchange/diffusion reactions
    jnz_before_filtering = np.where(sol_wt.fluxes.abs() > tolerance)[0]

    # Indices of exchange/diffusion reactions (not considered)
    eli_idx = [model.reactions.index(rxn_id) for rxn_id in elilist]

    jnz = np.setdiff1d(jnz_before_filtering, eli_idx)


    # Identify Single Lethal Reactions
    chunk_size = jnz.shape[0] // cpu_count()

    jsl_idx = list(
                   filter(
                          lambda rxn_idx: rxn_idx is not None,
                          Parallel(n_jobs=4,
                                   # threading performs better than
                                   # multiprocessing in only deletions
                                   backend='threading',
                                   verbose=5,
                                   batch_size=chunk_size)(
                                   delayed(_single_reactions_helper)
                                   (model, cutoff, gr_wt, i)
                                       for i in jnz)))

    # Indices -> Reactions
    jsl = model.reactions.get_by_any(jsl_idx)

    return jsl


def parallel_double_reactions(model, elilist, cutoff=0.01, **kwargs):
    '''
    Analyze double lethal reactions.

    Parameters
    ----------
    model: cobra.Model
        The model to perform operation on.
    elilist: list
        The list to be used for discarding the reactions.
    cutoff: float, optional (default 0.01)
        The value to be used as growth cutoff for deciding
        synthetic lethality.
    **kwargs

    Returns
    -------
    tuple(list of single synthetic lethal reaction IDs,
          list of double synthetic lethal reaction IDs)

    '''
    solver = 'glpk'
    tolerance = model.solver.configuration.tolerances.optimality

    model.solver = solver
    sol_wt = model.optimize()
    gr_wt = sol_wt.objective_value

    # Indices of non-zero flux reactions including exchange/diffusion reactions
    jnz_before_filtering = np.where(sol_wt.fluxes.abs() > tolerance)[0]

    # Indices of exchange/diffusion reactions (not considered)
    eli_idx = [model.reactions.index(rxn_id) for rxn_id in elilist]

    jnz = np.setdiff1d(jnz_before_filtering, eli_idx)

    # Identify single lethal reactions
    jsl = parallel_single_reactions(model, elilist, cutoff, solver=solver)

    # Indices of single lethal reacions
    jsl_idx = [model.reactions.index(jsl_id) for jsl_id in jsl]

    jnz_minus_jsl = np.setdiff1d(jnz, jsl_idx)

    jnz_minus_jsl_phase_2 = [rxn_pair for rxn_pair in product(jnz_minus_jsl,
                                                              repeat=2)]

    # Identify Double Lethal Reactions

    # Phase 1
    chunk_size_phase_1 = jnz_minus_jsl.shape[0] // cpu_count()
    jdl_idx_1 = list(
                     filter(
                            lambda rxn_pair_idx: rxn_pair_idx is not None,
                            Parallel(n_jobs=4,
                                     backend='multiprocessing',
                                     verbose=5,
                                     batch_size=chunk_size_phase_1)(
                                     delayed(_double_reactions_phase_1)
                                     (model, cutoff, tolerance, gr_wt, jnz,
                                      eli_idx, i) for i in jnz_minus_jsl)))

    # Phase 2
    chunk_size_phase_2 = len(jnz_minus_jsl_phase_2) // cpu_count()
    jdl_idx_2 = list(
                     filter(
                            lambda rxn_idx: rxn_idx is not None,
                            Parallel(n_jobs=4,
                                     backend='threading',
                                     verbose=5,
                                     batch_size=chunk_size_phase_2)(
                                     delayed(_double_reactions_phase_2)
                                     (model, cutoff, gr_wt, jnz_minus_jsl,
                                      pair)
                                         for pair in jnz_minus_jsl_phase_2)))

    jdl_idx = [inner_list for outer_list in list(filter(None, jdl_idx_1))
               for inner_list in outer_list] + jdl_idx_2

    # Indices -> Reactions
    jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in jdl_idx]

    return (jsl, jdl)
