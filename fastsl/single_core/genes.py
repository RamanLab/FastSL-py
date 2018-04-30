# -*- coding: utf-8 -*-

import math
from itertools import product

import numpy as np
from tqdm import tqdm


def single_sl_genes(model, cutoff, solver):
    '''Analysis for single lethal genes.'''

    model.solver = solver  # Basic solver configuration
    sol_wt = model.optimize()  # Identify minNorm flux distribution
    gr_wt = sol_wt.objective_value

    solver_tol = model.solver.configuration.tolerances.optimality

    # gives the non-zero flux reaction indices
    jnz_rxns = np.where(sol_wt.fluxes.abs() > solver_tol)[0]
    # set of the frozensets of genes associated with the jnz reactions
    # removes duplicates since set
    jnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                          for rxn_idx in jnz_rxns
                          if model.reactions[rxn_idx].gene_reaction_rule != ''}
    # indices of the genes obtained; removes duplicates too
    jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             jnz_rxns_genes_set
                                             for genes in single_frozenset]))

    # Identify Single Lethal Gene Deletions

    jsl_genes_idx = []

    for del_idx_i in tqdm(jnz_rxns_genes_idx, desc="Identifying jsl genes"):
        with model:
            model.genes[del_idx_i].knock_out()
            sol_ko_i = model.slim_optimize()
            if sol_ko_i < cutoff * gr_wt or math.isnan(sol_ko_i) is True:
                jsl_genes_idx.append(int(del_idx_i))

    # Indices -> Genes
    jsl_genes = model.genes.get_by_any(jsl_genes_idx)

    return jsl_genes


def double_sl_genes(model, cutoff, solver):
    '''Analysis for double lethal genes.'''

    model.solver = solver  # Basic solver configuration
    sol_wt = model.optimize()  # Identify minNorm flux distribution
    gr_wt = sol_wt.objective_value

    solver_tol = model.solver.configuration.tolerances.optimality

    # gives the non-zero flux reaction indices
    jnz_rxns = np.where(sol_wt.fluxes.abs() > solver_tol)[0]
    # set of the frozensets of genes associated with the jnz reactions
    # removes duplicates and non-gene associated reactions
    jnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                          for rxn_idx in jnz_rxns
                          if model.reactions[rxn_idx].gene_reaction_rule != ''}
    # indices of the genes obtained; removes duplicates too
    jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             jnz_rxns_genes_set
                                             for genes in single_frozenset]))

    # Identify single lethal genes
    jsl_genes = single_sl_genes(model,
                                cutoff,
                                solver)

    # Indices of single lethal genes
    jsl_genes_idx = [model.genes.index(jsl_id) for jsl_id in jsl_genes]

    jnz_copy = np.setdiff1d(jnz_rxns_genes_idx, jsl_genes_idx)  # jnz-jsl

    # Makes rxn pairs to remove nested for-loops in phase 2
    jnz_copy_phase_2 = [gene_pair for gene_pair in product(jnz_copy, repeat=2)]

    # Identify double lethal reactions

    jdl_genes_idx = []

    # Phase 1:
    for del_idx_i in tqdm(jnz_copy, desc="Identifying jdl genes: 1 of 2"):
        with model:
            model.genes[del_idx_i].knock_out()
            sol_ko_i = model.optimize()
            newnnz = np.where(sol_ko_i.fluxes.abs() > solver_tol)[0]
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

            for del_idx_j in jnz_i:
                with model:
                    model.genes[del_idx_j].knock_out()
                    sol_ko_ij = model.slim_optimize()
                    if sol_ko_ij < cutoff * gr_wt or \
                            math.isnan(sol_ko_ij) is True:
                        jdl_genes_idx.append([int(del_idx_i),
                                              int(del_idx_j)])

    # Phase 2:
    for del_idx_pair in tqdm(jnz_copy_phase_2,
                             desc="Identifying jdl genes: 2 of 2"):
        del_idx_i, del_idx_j = del_idx_pair
        if np.where(jnz_copy == del_idx_j) < np.where(jnz_copy == del_idx_i):
            with model:
                model.genes[del_idx_i].knock_out()
                model.genes[del_idx_j].knock_out()
                sol_ko_ij = model.slim_optimize()
                if sol_ko_ij < cutoff * gr_wt or \
                        math.isnan(sol_ko_ij) is True:
                    jdl_genes_idx.append([int(del_idx_i), int(del_idx_j)])

    # Indices -> Genes
    jdl_genes = [model.genes.get_by_any(gene_pair_idx) for gene_pair_idx
                 in jdl_genes_idx]

    return (jsl_genes, jdl_genes)
