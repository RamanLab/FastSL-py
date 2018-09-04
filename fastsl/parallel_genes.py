# -*- coding: utf-8 -*-

'''Contains functions to identify synthetic lethal genes.'''

import math
from os import cpu_count
from itertools import product

import numpy as np
from joblib import Parallel, delayed


def _single_genes_helper(model, cutoff, gr_wt, i):
    with model:
        model.genes[i].knock_out()
        sol_ko_i = model.slim_optimize()
        if sol_ko_i < cutoff * gr_wt or math.isnan(sol_ko_i) is True:
            return int(i)


def _double_genes_phase_1_helper(model, cutoff, gr_wt, i, j):
    with model:
        model.genes[i].knock_out()
        model.genes[j].knock_out()
        sol_ko_ij = model.slim_optimize()
        if sol_ko_ij < cutoff * gr_wt or \
                math.isnan(sol_ko_ij) is True:
            return [int(i), int(j)]


def _double_genes_phase_1(model, cutoff, gr_wt, tolerance, jnz, i):
    with model:
        model.genes[i].knock_out()
        sol_ko_i = model.optimize()
        newnnz = np.where(sol_ko_i.fluxes.abs() > tolerance)[0]
        newnnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                                 for rxn_idx in newnnz
                                 if model.reactions[rxn_idx]
                                 .gene_reaction_rule != ''}
        newnnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                          for single_frozenset in
                                          newnnz_rxns_genes_set
                                          for genes in single_frozenset]))

        jnz_i = np.setdiff1d(newnnz_rxns_genes_idx, jnz)

    return list(
                filter(
                       lambda rxn_pair_idx: rxn_pair_idx is not None,
                       Parallel(n_jobs=1,
                                backend='multiprocessing',
                                verbose=0,
                                batch_size='auto')(
                                delayed(_double_genes_phase_1_helper)
                                (model, cutoff, gr_wt, i, j)
                                    for j in jnz_i)))


def _double_lethal_genes_phase_2(model, cutoff, gr_wt, jnz_minus_jsl, pair):
    i, j = pair
    if np.where(jnz_minus_jsl == j) < np.where(jnz_minus_jsl == i):
        with model:
            model.genes[i].knock_out()
            model.genes[j].knock_out()
            sol_ko_ij = model.slim_optimize()
            if sol_ko_ij < cutoff * gr_wt or \
                    math.isnan(sol_ko_ij) is True:
                return [int(i), int(j)]


def parallel_single_genes(model, cutoff=0.01, **kwargs):
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
    tolerance = model.solver.configuration.tolerances.optimality

    model.solver = solver
    sol_wt = model.optimize()
    gr_wt = sol_wt.objective_value


    # Indices of the non-zero flux reactions
    jnz_rxns = np.where(sol_wt.fluxes.abs() > tolerance)[0]

    # Set of the frozensets of genes associated with the non-zero flux
    # reactions; removes duplicates
    jnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                          for rxn_idx in jnz_rxns
                          if model.reactions[rxn_idx].gene_reaction_rule != ''}

    # Indices of the genes obtained; removes duplicates
    jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             jnz_rxns_genes_set
                                             for genes in single_frozenset]))

    # Identify Single Lethal Genes
    chunk_size = jnz_rxns_genes_idx.shape[0] // cpu_count()

    jsl_genes_idx = list(
                         filter(
                                lambda rxn_idx: rxn_idx is not None,
                                Parallel(n_jobs=4,
                                         # threading performs better than
                                         # multiprocessing in only deletions
                                         backend='threading',
                                         verbose=5,
                                         batch_size=chunk_size)(
                                         delayed(_single_genes_helper)
                                         (model, cutoff, gr_wt, i)
                                             for i in jnz_rxns_genes_idx)))

    # Indices -> Genes
    jsl_genes = model.genes.get_by_any(jsl_genes_idx)

    return jsl_genes


def parallel_double_genes(model, cutoff=0.01, **kwargs):
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
    tolerance = model.solver.configuration.tolerances.optimality

    model.solver = solver
    sol_wt = model.optimize()
    gr_wt = sol_wt.objective_value


    # Indices of the non-zero flux reactions
    jnz_rxns = np.where(sol_wt.fluxes.abs() > tolerance)[0]

    # Set of the frozensets of genes associated with the non-zero flux
    # reactions; removes duplicates
    jnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                          for rxn_idx in jnz_rxns
                          if model.reactions[rxn_idx].gene_reaction_rule != ''}

    # Indices of the genes obtained; removes duplicates
    jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             jnz_rxns_genes_set
                                             for genes in single_frozenset]))

    # Identify single lethal genes
    jsl_genes = parallel_single_genes(model, cutoff, solver=solver)

    # Indices of single lethal genes
    jsl_genes_idx = [model.genes.index(jsl_id) for jsl_id in jsl_genes]

    jnz_minus_jsl = np.setdiff1d(jnz_rxns_genes_idx, jsl_genes_idx)

    jnz_minus_jsl_phase_2 = [rxn_pair for rxn_pair in product(jnz_minus_jsl,
                                                         repeat=2)]

    # Identify Double Lethal Genes

    # Phase 1
    chunk_size_phase_1 = jnz_minus_jsl.shape[0] // cpu_count()
    jdl_idx_1 = list(
                     filter(
                            lambda rxn_pair_idx: rxn_pair_idx is not None,
                            Parallel(n_jobs=4,
                                     backend='multiprocessing',
                                     verbose=5,
                                     batch_size=chunk_size_phase_1)(
                                     delayed(_double_genes_phase_1)
                                     (model, cutoff, gr_wt, tolerance,
                                      jnz_rxns_genes_idx, i)
                                         for i in jnz_minus_jsl)))

    # Phase 2
    chunk_size_phase_2 = len(jnz_minus_jsl_phase_2) // cpu_count()
    jdl_idx_2 = list(
                     filter(
                            lambda rxn_idx: rxn_idx is not None,
                            Parallel(n_jobs=4,
                                     backend='threading',
                                     verbose=5,
                                     batch_size=chunk_size_phase_2)(
                                     delayed(_double_lethal_genes_phase_2)
                                     (model, cutoff, gr_wt, jnz_minus_jsl,
                                      pair)
                                         for pair in jnz_minus_jsl_phase_2)))

    jdl_idx = [inner_list for outer_list in list(filter(None, jdl_idx_1))
               for inner_list in outer_list] + jdl_idx_2

    # Indices -> Genes
    jdl_genes = [model.genes.get_by_any(rxn_pair_idx) for rxn_pair_idx
                 in jdl_idx]

    return (jsl_genes, jdl_genes)
