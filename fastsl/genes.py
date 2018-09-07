# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

'''Contains functions to identify double synthetic lethal deletions.'''

import math
from itertools import product

import numpy as np
from joblib import Parallel, delayed


def _single_genes_helper(model, tolerance, idx):
    with model:
        model.genes[idx].knock_out()
        sol_ko_i = model.slim_optimize()
        if sol_ko_i < tolerance or math.isnan(sol_ko_i) is True:
            return int(idx)


def single_genes(model, cutoff=0.01, **kwargs):
    '''
    Analyze single lethal genes.

    Parameters
    ----------
    model = cobra.Model
        The model to perform operation on.
    cutoff: float, optional (default 0.01)
        The value to be used as growth cutoff for deciding
        synthetic lethality.
    **kwargs

    Returns
    -------
    list of single lethal gene IDs

    '''
    solver = 'glpk'

    model.solver = solver
    sol = model.optimize()
    obj = sol.objective_value

    tolerance = cutoff * obj

    # Obtain indices of the non-zero flux reactions (including
    # exchange/diffusion reactions)
    jnz_rxns = np.where(sol.fluxes.abs() > tolerance)[0]

    # Obtain set of the frozensets of genes associated with the non-zero flux
    # reactions; removes duplicates
    jnz_rxns_genes = {model.reactions[rxn_idx].genes
                      for rxn_idx in jnz_rxns
                      if model.reactions[rxn_idx].gene_reaction_rule != ''}

    # Obtain indices of the genes obtained; removes duplicates
    genes_idx = np.unique(np.array([model.genes.index(genes)
                                    for single_frozenset in
                                    jnz_rxns_genes
                                    for genes in single_frozenset]))

    # Identify Single Lethal Genes

    jsl_genes_idx = list(filter(lambda rxn_idx: rxn_idx is not None,
                                Parallel(n_jobs=-1,
                                         backend='threading',
                                         verbose=1)(
                                             delayed(_single_genes_helper)
                                             (model, tolerance, idx)
                                             for idx in genes_idx)))


    # Indices -> Genes
    jsl_genes = model.genes.get_by_any(jsl_genes_idx)

    return jsl_genes


def _double_genes_phase_1_helper(model, tolerance, i_idx, j_idx):
    with model:
        model.genes[i_idx].knock_out()
        model.genes[j_idx].knock_out()
        sol_ko_ij = model.slim_optimize()
        if sol_ko_ij < tolerance or math.isnan(sol_ko_ij) is True:
            return [int(i_idx), int(j_idx)]


def _double_genes_phase_1(model, tolerance, jnz, i_idx):
    with model:
        model.genes[i_idx].knock_out()
        sol_ko_i = model.optimize()
        new_nz = np.where(sol_ko_i.fluxes.abs() > 1E-6)[0]
        new_nz_rxns_genes = {model.reactions[rxn_idx].genes
                             for rxn_idx in new_nz
                             if model.reactions[rxn_idx]
                             .gene_reaction_rule != ''}
        new_nz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                    for single_frozenset in
                                                    new_nz_rxns_genes
                                                    for genes in single_frozenset]))

        jnz_i = np.setdiff1d(new_nz_rxns_genes_idx, jnz)

    return list(filter(lambda rxn_pair_idx: rxn_pair_idx is not None,
                       Parallel(n_jobs=-1,
                                backend='threading',
                                verbose=1)(
                                    delayed(_double_genes_phase_1_helper)
                                    (model, tolerance, i_idx, j_idx)
                                    for j_idx in jnz_i)))


def _double_lethal_genes_phase_2(model, tolerance, jnz_minus_jsl, pair):
    i_idx, j_idx = pair
    if np.where(jnz_minus_jsl == j_idx) < np.where(jnz_minus_jsl == i_idx):
        with model:
            model.genes[i_idx].knock_out()
            model.genes[j_idx].knock_out()
            sol_ko_ij = model.slim_optimize()
            if sol_ko_ij < tolerance or math.isnan(sol_ko_ij) is True:
                return [int(i_idx), int(j_idx)]


def double_genes(model, cutoff=0.01, **kwargs):
    '''
    Analyze double lethal genes.

    Parameters
    ----------
    model: cobra.Model
        The model to perform operation on.
    cutoff: float, optional (default 0.01)
        The value to be used as growth cutoff for deciding
        synthetic lethality.
    **kwargs

    Returns
    -------
    tuple(list of single lethal gene IDs,
          list of double lethal gene IDs)

    '''
    solver = 'glpk'

    model.solver = solver
    sol = model.optimize()
    obj = sol.objective_value


    tolerance = cutoff * obj

    # Obtain indices of the non-zero flux reactions (including
    # exchange/diffusion reactions)
    jnz_rxns = np.where(sol.fluxes.abs() > tolerance)[0]

    # Obtain set of the frozensets of genes associated with the non-zero flux
    # reactions; removes duplicates
    jnz_rxns_genes = {model.reactions[rxn_idx].genes
                      for rxn_idx in jnz_rxns
                      if model.reactions[rxn_idx].gene_reaction_rule != ''}

    # Obtain indices of the genes obtained; removes duplicates
    jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             jnz_rxns_genes
                                             for genes in single_frozenset]))

    # Identify Single Lethal Genes
    jsl_genes = single_genes(model, cutoff, solver=solver)

    # Indices of Single Lethal Genes
    jsl_genes_idx = [model.genes.index(jsl_id) for jsl_id in jsl_genes]

    jnz_minus_jsl = np.setdiff1d(jnz_rxns_genes_idx, jsl_genes_idx)

    # Make reaction pairs to remove nested for-loops in phase 2
    jnz_minus_jsl_phase_2 = [gene_pair for gene_pair in product(jnz_minus_jsl,
                                                           repeat=2)]

    # Identify Double Lethal Genes
    
    # Phase 1:
    jdl_idx_1 = list(filter(lambda rxn_pair_idx: rxn_pair_idx is not None,
                            Parallel(n_jobs=-1,
                                     backend='threading',
                                     verbose=1)(
                                         delayed(_double_genes_phase_1)
                                         (model, tolerance,
                                          jnz_rxns_genes_idx, idx)
                                         for idx in jnz_minus_jsl)))


    # for i in tqdm(jnz_minus_jsl, desc=("Identifying double lethal genes:"
    #                                    " 1 of 2")):
    #     with model:
    #         model.genes[i].knock_out()
    #         sol_ko_i = model.optimize()
    #         newnnz = np.where(sol_ko_i.fluxes.abs() > tolerance)[0]
    #         newnnz_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
    #                                        for rxn_idx in newnnz
    #                                        if model.reactions[rxn_idx]
    #                                        .gene_reaction_rule != ''}
    #         newnnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
    #                                                     for single_frozenset
    #                                                     in newnnz_rxns_genes_frozenset
    #                                                     for genes in single_frozenset]))
    #         jnz_i = np.setdiff1d(newnnz_rxns_genes_idx,
    #                              jnz_rxns_genes_idx)

    #         for j in jnz_i:
    #             with model:
    #                 model.genes[j].knock_out()
    #                 sol_ko_ij = model.slim_optimize()
    #                 if sol_ko_ij < cutoff * gr_wt or \
    #                         math.isnan(sol_ko_ij) is True:
    #                     jdl_genes_idx.append([int(i),
    #                                           int(j)])

    # Phase 2:
    jdl_idx_2 = list(filter(lambda rxn_idx: rxn_idx is not None,
                            Parallel(n_jobs=-1,
                                     backend='threading',
                                     verbose=1)(
                                         delayed(_double_lethal_genes_phase_2)
                                         (model, tolerance, jnz_minus_jsl,
                                          pair)
                                         for pair in jnz_minus_jsl_phase_2)))


    # for pair in tqdm(jnz_minus_jsl_phase_2,
    #                  desc=("Identifying double lethal genes:"
    #                        " 2 of 2")):
    #     i, j = pair
    #     if np.where(jnz_minus_jsl == j) < np.where(jnz_minus_jsl == i):
    #         with model:
    #             model.genes[i].knock_out()
    #             model.genes[j].knock_out()
    #             sol_ko_ij = model.slim_optimize()
    #             if sol_ko_ij < cutoff * gr_wt or \
    #                     math.isnan(sol_ko_ij) is True:
    #                 jdl_genes_idx.append([int(i), int(j)])

    jdl_genes_idx = [inner_list for outer_list in list(filter(None, jdl_idx_1))
               for inner_list in outer_list] + jdl_idx_2

    # Indices -> Genes
    jdl_genes = [model.genes.get_by_any(gene_pair_idx) for gene_pair_idx
                 in jdl_genes_idx]

    return (jsl_genes, jdl_genes)
