# -*- coding: utf-8 -*-

import math
from itertools import product

import numpy as np
from tqdm import tqdm


def single_sl(model, cutoff, eli_list, solver):
    '''Analysis for single lethal reactions.'''
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

    jsl_idx = []

    for del_idx_i in tqdm(jnz, desc="Identifying jsl reactions"):
        with model:
            model.reactions[del_idx_i].knock_out()
            sol_ko_i = model.slim_optimize()
            if sol_ko_i < cutoff * gr_wt or math.isnan(sol_ko_i) is True:
                jsl_idx.append(int(del_idx_i))

    # Indices -> Reactions
    jsl = model.reactions.get_by_any(jsl_idx)
    # jsl = [model.reactions[rxn_idx].id for rxn_idx in jsl_idx]

    return jsl


def double_sl(model, cutoff, eli_list, solver):
    '''Analysis for double lethal reactions.'''
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
    jsl = single_sl(model,
                    cutoff,
                    eli_list,
                    solver)

    # Indices of single lethal reacions
    jsl_idx = [model.reactions.index(jsl_id) for jsl_id in jsl]

    jnz_copy = np.setdiff1d(jnz, jsl_idx)  # jnz-jsl

    # Makes rxn pairs to remove nested for-loops in phase 2
    jnz_copy_phase_2 = [rxn_pair for rxn_pair in product(jnz_copy, repeat=2)]

    # Identify double lethal reactions

    jdl_idx = []

    # Phase 1:
    for del_idx_i in tqdm(jnz_copy, desc="Identifying jdl reactions: 1 of 2"):
        with model:
            model.reactions[del_idx_i].knock_out()
            sol_ko_i = model.optimize()
            newnnz = np.where(sol_ko_i.fluxes.abs() > solver_tol)[0]
            jnz_i_before_filtering = np.setdiff1d(newnnz, jnz)
            jnz_i = np.setdiff1d(jnz_i_before_filtering, eli_idx)

            for del_idx_j in jnz_i:
                with model:
                    model.reactions[del_idx_j].knock_out()
                    sol_ko_ij = model.slim_optimize()
                    if sol_ko_ij < cutoff * gr_wt or \
                            math.isnan(sol_ko_ij) is True:
                        jdl_idx.append([int(del_idx_i), int(del_idx_j)])

    # Phase 2:
    for del_idx_pair in tqdm(jnz_copy_phase_2,
                             desc="Identifying jdl reactions: 2 of 2"):
        del_idx_i, del_idx_j = del_idx_pair
        if np.where(jnz_copy == del_idx_j) < np.where(jnz_copy == del_idx_i):
            with model:
                model.reactions[del_idx_i].knock_out()
                model.reactions[del_idx_j].knock_out()
                sol_ko_ij = model.slim_optimize()
                if sol_ko_ij < cutoff * gr_wt or \
                        math.isnan(sol_ko_ij) is True:
                    jdl_idx.append([int(del_idx_i), int(del_idx_j)])

    # Indices -> Reactions
    jdl = [model.reactions.get_by_any(rxn_pair_idx) for rxn_pair_idx
           in jdl_idx]

    return (jsl, jdl)
