# -*- coding: utf-8 -*-

'''Contains functions to identify single synthetic lethal reactions.'''

import math
from itertools import product

import numpy as np

from tqdm import tqdm


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
    sol_wt = model.optimize()
    gr_wt = sol_wt.objective_value


    # Indices of non-zero flux reactions including exchange/diffusion reactions
    jnz_before_filtering = np.where(sol_wt.fluxes.abs() > tolerance)[0]

    # Indices of exchange/diffusion reactions (not considered)
    eli_idx = [model.reactions.index(rxn_id) for rxn_id in elilist]

    jnz = np.setdiff1d(jnz_before_filtering, eli_idx)

    # Identify Single Lethal Reaction Deletions
    jsl_idx = []

    for i in tqdm(jnz, desc="Identifying single lethal reactions:"):
        with model:
            model.reactions[i].knock_out()
            sol_ko_i = model.slim_optimize()
            if sol_ko_i < cutoff * gr_wt or math.isnan(sol_ko_i) is True:
                jsl_idx.append(int(i))

    # Indices -> Reactions
    jsl = model.reactions.get_by_any(jsl_idx)

    return jsl


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
    jsl = single_reactions(model, elilist, cutoff, solver=solver)

    # Indices of single lethal reacions
    jsl_idx = [model.reactions.index(jsl_id) for jsl_id in jsl]

    jnz_minus_jsl = np.setdiff1d(jnz, jsl_idx)

    # Make reaction pairs to remove nested for-loops in phase 2
    jnz_minus_jsl_phase_2 = [rxn_pair for rxn_pair in product(jnz_minus_jsl,
                                                              repeat=2)]

    # Identify double lethal reactions
    jdl_idx = []

    # Phase 1:
    for i in tqdm(jnz_minus_jsl, desc=("Identifying double lethal reactions:"
                                       " 1 of 2")):
        with model:
            model.reactions[i].knock_out()
            sol_ko_i = model.optimize()
            newnnz = np.where(sol_ko_i.fluxes.abs() > tolerance)[0]
            jnz_i_before_filtering = np.setdiff1d(newnnz, jnz)
            jnz_i = np.setdiff1d(jnz_i_before_filtering, eli_idx)

            for j in jnz_i:
                with model:
                    model.reactions[j].knock_out()
                    sol_ko_ij = model.slim_optimize()
                    if sol_ko_ij < cutoff * gr_wt or \
                            math.isnan(sol_ko_ij) is True:
                        jdl_idx.append([int(i), int(j)])

    # Phase 2:
    for pair in tqdm(jnz_minus_jsl_phase_2,
                     desc=("Identifying double lethal reactions:"
                           " 2 of 2")):
        i, j = pair
        if np.where(jnz_minus_jsl == j) < np.where(jnz_minus_jsl == i):
            with model:
                model.reactions[i].knock_out()
                model.reactions[j].knock_out()
                sol_ko_ij = model.slim_optimize()
                if sol_ko_ij < cutoff * gr_wt or \
                        math.isnan(sol_ko_ij) is True:
                    jdl_idx.append([int(i), int(j)])

    # Indices -> Reactions
    jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in jdl_idx]

    return (jsl, jdl)
