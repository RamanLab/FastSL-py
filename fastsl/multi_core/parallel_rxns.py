# -*- coding: utf-8 -*-

import math
from os import cpu_count
from itertools import product

import numpy as np
from joblib import Parallel, delayed


def _single_lethal_reactions(model, cutoff, gr_wt, del_idx_i):
    with model:
        model.reactions[del_idx_i].knock_out()
        sol_ko_i = model.slim_optimize()
        if sol_ko_i < cutoff * gr_wt or math.isnan(sol_ko_i) is True:
            return int(del_idx_i)
        else:
            return None


def _double_lethal_reactions_phase_1_helper(model, cutoff, gr_wt, del_idx_i,
                                            del_idx_j):
    with model:
        model.reactions[del_idx_i].knock_out()
        model.reactions[del_idx_j].knock_out()
        sol_ko_ij = model.slim_optimize()
        if sol_ko_ij < cutoff * gr_wt or \
                math.isnan(sol_ko_ij) is True:
            return [int(del_idx_i), int(del_idx_j)]
        else:
            return None


def _double_lethal_reactions_phase_1(model, solver_tol, cutoff, gr_wt,
                                     jnz, eli_idx, del_idx_i):
    with model:
        model.reactions[del_idx_i].knock_out()
        sol_ko_i = model.optimize()
        newnnz = np.where(sol_ko_i.fluxes.abs() > solver_tol)[0]
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
                                    _double_lethal_reactions_phase_1_helper)
                                (model,
                                 cutoff,
                                 gr_wt,
                                 del_idx_i,
                                 del_idx_j) for del_idx_j in jnz_i)))


def _double_lethal_reactions_phase_2(model, cutoff, gr_wt, jnz_copy,
                                     del_idx_pair):
    del_idx_i, del_idx_j = del_idx_pair
    if np.where(jnz_copy == del_idx_j) < np.where(jnz_copy == del_idx_i):
        with model:
            model.reactions[del_idx_i].knock_out()
            model.reactions[del_idx_j].knock_out()
            sol_ko_ij = model.slim_optimize()
            if sol_ko_ij < cutoff * gr_wt or \
                    math.isnan(sol_ko_ij) is True:
                return [int(del_idx_i), int(del_idx_j)]
            else:
                return None


def parallel_single_sl(model, cutoff, eli_list, solver):
    ''' Analysis for single lethal reactions.'''
    model.solver = solver  # Basic solver configuration
    sol_wt = model.optimize()  # Identify minNorm flux distribution
    gr_wt = sol_wt.objective_value

    solver_tol = model.solver.configuration.tolerances.optimality

    # Indices of non-zero flux reactions including exchange/diffusion reactions
    jnz_before_filtering = np.where(sol_wt.fluxes.abs() > solver_tol)[0]

    # Indices of exchange/diffusion reactions (not considered)
    eli_idx = [model.reactions.index(reaction_id) for reaction_id in eli_list]

    jnz = np.setdiff1d(jnz_before_filtering, eli_idx)  # jnz

    # Identify Single Lethal Reaction Deletions

    chunk_size = jnz.shape[0] // cpu_count()  # integer division

    jsl_idx = list(
                   filter(
                          lambda rxn_idx: rxn_idx is not None,
                          Parallel(n_jobs=4,
                                   # threading performs better than
                                   # multiprocessing in only deletions
                                   backend='threading',
                                   verbose=5,
                                   batch_size=chunk_size)(
                                   delayed(_single_lethal_reactions)
                                   (model,
                                    cutoff,
                                    gr_wt,
                                    del_idx_i) for del_idx_i in jnz)))

    # Indices -> Reactions
    jsl = model.reactions.get_by_any(jsl_idx)

    return jsl


def parallel_double_sl(model, cutoff, eli_list, solver):
    ''' Analysis for double lethal reactions '''
    model.solver = solver  # Basic solver configuration
    sol_wt = model.optimize()  # Identify minNorm flux distribution
    gr_wt = sol_wt.objective_value

    solver_tol = model.solver.configuration.tolerances.optimality

    # Indices of non-zero flux reactions including exchange/diffusion reactions
    jnz_before_filtering = np.where(sol_wt.fluxes.abs() > solver_tol)[0]

    # Indices of exchange/diffusion reactions (not considered)
    eli_idx = [model.reactions.index(reaction_id) for reaction_id in eli_list]

    jnz = np.setdiff1d(jnz_before_filtering, eli_idx)  # jnz

    # Identify single lethal reactions
    jsl = parallel_single_sl(model,
                             cutoff,
                             eli_list,
                             solver)

    # Indices of single lethal reacions
    jsl_idx = [model.reactions.index(jsl_id) for jsl_id in jsl]

    jnz_copy = np.setdiff1d(jnz, jsl_idx)  # jnz-jsl

    jnz_copy_phase_2 = [rxn_pair for rxn_pair in product(jnz_copy, repeat=2)]

    # Identify Double Lethal Reaction Deletions

    # Phase 1
    chunk_size_phase_1 = jnz_copy.shape[0] // cpu_count()
    jdl_idx_1 = list(
                     filter(
                            lambda rxn_pair_idx: rxn_pair_idx is not None,
                            Parallel(n_jobs=4,
                                     backend='multiprocessing',
                                     verbose=5,
                                     batch_size=chunk_size_phase_1)(
                                     delayed(_double_lethal_reactions_phase_1)
                                     (model,
                                      solver_tol,
                                      cutoff,
                                      gr_wt,
                                      jnz,
                                      eli_idx,
                                      del_idx_i) for del_idx_i in jnz_copy)))

    # Phase 2
    chunk_size_phase_2 = len(jnz_copy_phase_2) // cpu_count()
    jdl_idx_2 = list(
                     filter(
                            lambda rxn_idx: rxn_idx is not None,
                            Parallel(n_jobs=4,
                                     backend='threading',
                                     verbose=5,
                                     batch_size=chunk_size_phase_2)(
                                     delayed(_double_lethal_reactions_phase_2)
                                     (model,
                                      cutoff,
                                      gr_wt,
                                      jnz_copy,
                                      del_idx_pair) for del_idx_pair in jnz_copy_phase_2)))

    jdl_idx = [inner_list for outer_list in list(filter(None, jdl_idx_1))
               for inner_list in outer_list] + jdl_idx_2

    # Indices -> Reactions
    jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in jdl_idx]

    return (jsl, jdl)
