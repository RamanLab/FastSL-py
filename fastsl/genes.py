# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

'''Contains functions to identify double synthetic lethal deletions.'''

import math
from itertools import product

import numpy as np

from tqdm import tqdm


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
    jsl_genes_idx = []

    for i in tqdm(jnz_rxns_genes_idx,
                  desc="Identifying single lethal genes: "):
        with model:
            model.genes[i].knock_out()
            sol_ko_i = model.slim_optimize()
            if sol_ko_i < cutoff * gr_wt or math.isnan(sol_ko_i) is True:
                jsl_genes_idx.append(int(i))

    # Indices -> Genes
    jsl_genes = model.genes.get_by_any(jsl_genes_idx)

    return jsl_genes


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
    jsl_genes = single_genes(model, cutoff, solver=solver)

    # Indices of single lethal genes
    jsl_genes_idx = [model.genes.index(jsl_id) for jsl_id in jsl_genes]

    jnz_minus_jsl = np.setdiff1d(jnz_rxns_genes_idx, jsl_genes_idx)

    # Make reaction pairs to remove nested for-loops in phase 2
    jnz_minus_jsl_phase_2 = [gene_pair for gene_pair in product(jnz_minus_jsl,
                                                           repeat=2)]

    # Identify Double Lethal Genes
    jdl_genes_idx = []

    # Phase 1:
    for i in tqdm(jnz_minus_jsl, desc=("Identifying double lethal genes:"
                                       " 1 of 2")):
        with model:
            model.genes[i].knock_out()
            sol_ko_i = model.optimize()
            newnnz = np.where(sol_ko_i.fluxes.abs() > tolerance)[0]
            newnnz_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                           for rxn_idx in newnnz
                                           if model.reactions[rxn_idx]
                                           .gene_reaction_rule != ''}
            newnnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                        for single_frozenset
                                                        in newnnz_rxns_genes_frozenset
                                                        for genes in single_frozenset]))
            jnz_i = np.setdiff1d(newnnz_rxns_genes_idx,
                                 jnz_rxns_genes_idx)

            for j in jnz_i:
                with model:
                    model.genes[j].knock_out()
                    sol_ko_ij = model.slim_optimize()
                    if sol_ko_ij < cutoff * gr_wt or \
                            math.isnan(sol_ko_ij) is True:
                        jdl_genes_idx.append([int(i),
                                              int(j)])

    # Phase 2:
    for pair in tqdm(jnz_minus_jsl_phase_2,
                     desc=("Identifying double lethal genes:"
                           " 2 of 2")):
        i, j = pair
        if np.where(jnz_minus_jsl == j) < np.where(jnz_minus_jsl == i):
            with model:
                model.genes[i].knock_out()
                model.genes[j].knock_out()
                sol_ko_ij = model.slim_optimize()
                if sol_ko_ij < cutoff * gr_wt or \
                        math.isnan(sol_ko_ij) is True:
                    jdl_genes_idx.append([int(i), int(j)])

    # Indices -> Genes
    jdl_genes = [model.genes.get_by_any(gene_pair_idx) for gene_pair_idx
                 in jdl_genes_idx]

    return (jsl_genes, jdl_genes)
