# -*- coding: utf-8 -*-

'''Contains functions to identify synthetic lethal reactions.'''

import math
from itertools import product

import numpy as np
from joblib import Parallel, delayed


def _single_reactions_helper(model, tolerance, idx):
    with model:
        model.reactions[idx].knock_out()
        sol_ko_i = model.slim_optimize()
        if sol_ko_i < tolerance or math.isnan(sol_ko_i) is True:
            return int(idx)


def single_reactions(model, elilist, cutoff=0.01, **kwargs):
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
    sol = model.optimize()
    obj = sol.objective_value

    tolerance = cutoff * obj

    # Obtain indices of non-zero flux reactions (including
    # exchange/diffusion reactions)
    jnz = np.where(sol.fluxes.abs() > 1E-6)[0]

    # Remove indices of exchange/diffusion reactions
    eli_idx = [model.reactions.index(rxn_id) for rxn_id in elilist]
    jnz = np.setdiff1d(jnz, eli_idx)

    # Identify Single Lethal Reactions
    jsl_idx = list(filter(lambda idx: idx is not None,
                          Parallel(n_jobs=-1,
                                   backend='threading',
                                   verbose=1)(
                                       delayed(_single_reactions_helper)
                                       (model, tolerance, idx)
                                       for idx in jnz)))

    # Indices -> Reactions
    jsl = model.reactions.get_by_any(jsl_idx)

    return jsl


def _double_reactions_phase_1_helper(model, tolerance, i_idx, j_idx):
    with model:
        model.reactions[i_idx].knock_out()
        model.reactions[j_idx].knock_out()
        sol_ko_ij = model.slim_optimize()
        if sol_ko_ij < tolerance or math.isnan(sol_ko_ij) is True:
            return [int(i_idx), int(j_idx)]


def _double_reactions_phase_1(model, tolerance, jnz, eli_idx, idx):
    with model:
        model.reactions[idx].knock_out()
        sol_ko_i = model.optimize()
        new_nz = np.where(sol_ko_i.fluxes.abs() > 1E-6)[0]
        jnz_i = np.setdiff1d(new_nz, jnz)
        jnz_i = np.setdiff1d(jnz_i, eli_idx)

    return list(filter(lambda rxn_pair_idx: rxn_pair_idx is not None,
                       Parallel(n_jobs=-1,
                                backend='threading')(
                                    delayed(_double_reactions_phase_1_helper)
                                    (model, tolerance, idx, j_idx)
                                    for j_idx in jnz_i)))


def _double_reactions_phase_2(model, tolerance, jnz_minus_jsl, pair):
    i_idx, j_idx = pair
    if np.where(jnz_minus_jsl == j_idx) < np.where(jnz_minus_jsl == i_idx):
        with model:
            model.reactions[i_idx].knock_out()
            model.reactions[j_idx].knock_out()
            sol_ko_ij = model.slim_optimize()
            if sol_ko_ij < tolerance or math.isnan(sol_ko_ij) is True:
                return [int(i_idx), int(j_idx)]


def double_reactions(model, elilist, cutoff=0.01, **kwargs):
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

    model.solver = solver
    sol = model.optimize()
    obj = sol.objective_value

    tolerance = cutoff * obj

    # Obtain indices of non-zero flux reactions (including exchange/diffusion
    # reactions)
    jnz = np.where(sol.fluxes.abs() > 1E-6)[0]

    # Remove indices of exchange/diffusion reactions
    eli_idx = [model.reactions.index(rxn_id) for rxn_id in elilist]

    jnz = np.setdiff1d(jnz, eli_idx)

    # Identify Single Lethal Reactions
    jsl = single_reactions(model, elilist, cutoff, solver=solver)

    # Indices of Single Lethal Reactions
    jsl_idx = [model.reactions.index(jsl_id) for jsl_id in jsl]

    jnz_minus_jsl = np.setdiff1d(jnz, jsl_idx)

    # Make reaction pairs to remove nested for-loops in phase 2
    jnz_minus_jsl_phase_2 = [rxn_pair for rxn_pair in product(jnz_minus_jsl,
                                                              repeat=2)]

    # Identify Double Lethal Reactions
    # jdl_idx = []

    # Phase 1:
    jdl_idx_1 = list(filter(lambda rxn_pair_idx: rxn_pair_idx is not None,
                            Parallel(n_jobs=-1,
                                     backend='threading',
                                     verbose=1)(
                                         delayed(_double_reactions_phase_1)
                                         (model, tolerance, jnz, eli_idx, idx)
                                         for idx in jnz_minus_jsl)))

    # for i in tqdm(jnz_minus_jsl, desc=("Identifying double lethal reactions:"
    #                                    " 1 of 2")):
    #     with model:
    #         model.reactions[i].knock_out()
    #         sol_ko_i = model.optimize()
    #         newnnz = np.where(sol_ko_i.fluxes.abs() > tolerance)[0]
    #         jnz_i = np.setdiff1d(newnnz, jnz)
    #         jnz_i = np.setdiff1d(jnz_i, eli_idx)

    #         for j in jnz_i:
    #             with model:
    #                 model.reactions[j].knock_out()
    #                 sol_ko_ij = model.slim_optimize()
    #                 if sol_ko_ij < cutoff * gr_wt or \
    #                         math.isnan(sol_ko_ij) is True:
    #                     jdl_idx.append([int(i), int(j)])

    # Phase 2:

    jdl_idx_2 = list(filter(lambda rxn_idx: rxn_idx is not None,
                            Parallel(n_jobs=-1,
                                     backend='threading',
                                     verbose=1)(
                                         delayed(_double_reactions_phase_2)
                                         (model, tolerance, jnz_minus_jsl,
                                          pair)
                                         for pair in jnz_minus_jsl_phase_2)))

    # for pair in tqdm(jnz_minus_jsl_phase_2,
    #                  desc=("Identifying double lethal reactions:"
    #                        " 2 of 2")):
    #     i, j = pair
    #     if np.where(jnz_minus_jsl == j) < np.where(jnz_minus_jsl == i):
    #         with model:
    #             model.reactions[i].knock_out()
    #             model.reactions[j].knock_out()
    #             sol_ko_ij = model.slim_optimize()
    #             if sol_ko_ij < cutoff * gr_wt or \
    #                     math.isnan(sol_ko_ij) is True:
    #                 jdl_idx.append([int(i), int(j)])

    jdl_idx = [inner_list for outer_list in list(filter(None, jdl_idx_1))
               for inner_list in outer_list] + jdl_idx_2

    # Indices -> Reactions
    jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in jdl_idx]

    return (jsl, jdl)
