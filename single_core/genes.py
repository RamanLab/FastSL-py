# -*- coding: utf-8 -*-

import math
from itertools import product

import numpy as np
from tqdm import tqdm


def single_sl_genes(model, cutoff, solver):
    """
    Analysis for single lethal genes
    """

    model.solver = solver  # Basic solver configuration
    sol_wt = model.optimize()  # Identify minNorm flux distribution
    gr_wt = sol_wt.objective_value

    # gives the non-zero flux reaction indices
    jnz_rxns = np.flatnonzero(sol_wt.fluxes)
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
    """
    Analysis for double lethal genes
    """

    model.solver = solver  # Basic solver configuration
    sol_wt = model.optimize()  # Identify minNorm flux distribution
    gr_wt = sol_wt.objective_value

    # gives the non-zero flux reaction indices
    jnz_rxns = np.flatnonzero(sol_wt.fluxes)
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
            newnnz = np.flatnonzero(sol_ko_i.fluxes)
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


def triple_sl_genes(model, cutoff, solver):
    """
    Analysis for triple lethal reactions
    """

    model.solver = solver  # Basic solver configuration
    sol_wt = model.optimize()  # Identify minNorm flux distribution
    gr_wt = sol_wt.objective_value

    # gives the non-zero flux reaction indices
    jnz_rxns = np.flatnonzero(sol_wt.fluxes)
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

    # Identify triple lethal reactions

    jdl_genes_idx = []
    jtl_genes_idx = []

    # Phase 1:
    for del_idx_i in tqdm(jnz_copy,
                          desc="Identifying jdl and jtl genes: 1 of 2"):
        with model:
            model.genes[del_idx_i].knock_out()
            sol_ko_i = model.optimize()
            newnnz = np.flatnonzero(sol_ko_i.fluxes)
            newnnz_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                           for rxn_idx in newnnz
                                           if model.reactions[rxn_idx]
                                           .gene_reaction_rule != ''}
            newnnz_rxns_genes_idx = np.unique(np.array([model.genes.
                                                        index(genes)
                                                        for single_frozenset
                                                        in newnnz_rxns_genes_frozenset
                                                        for genes in single_frozenset]))
            jnz_i = np.setdiff1d(newnnz_rxns_genes_idx,
                                 jnz_rxns_genes_idx)

            for del_idx_j in jnz_i:
                with model:
                    model.genes[del_idx_j].knock_out()
                    sol_ko_ij = model.optimize()
                    if sol_ko_ij.objective_value < cutoff * gr_wt or \
                            sol_ko_ij.status != 'infeasible':
                        jdl_genes_idx.append([int(del_idx_i),
                                              int(del_idx_j)])

                    elif sol_ko_ij.status == 'infeasible':
                        sol_ko_ij = model.optimize()
                        if sol_ko_ij.objective_value < cutoff * gr_wt and \
                                math.isnan(sol_ko_ij.objective_value) is True:
                            jdl_genes_idx.append([int(del_idx_i),
                                                  int(del_idx_j)])
                            # continue

                        newnnz_ij = np.flatnonzero(sol_ko_ij.fluxes)
                        newnnz_ij_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                                          for rxn_idx in newnnz
                                                          if model.reactions[rxn_idx].gene_reaction_rule != ''}
                        newnnz_ij_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                                       for single_frozenset
                                                                       in newnnz_rxns_genes_frozenset
                                                                       for genes in single_frozenset]))
                        jnz_ij = np.setdiff1d(newnnz_ij_rxns_genes_idx,
                                              jnz_rxns_genes_idx)
                        for del_idx_k in jnz_ij:
                            with model:
                                model.genes[del_idx_k].knock_out()
                                sol_ko_ijk = model.slim_optimize()
                                if sol_ko_ijk < cutoff * gr_wt or \
                                        math.isnan(sol_ko_ijk) is True:
                                    jtl_genes_idx.append([int(del_idx_i),
                                                          int(del_idx_j),
                                                          int(del_idx_k)])

    # Phase 2:
    for del_idx_pair in tqdm(jnz_copy_phase_2,
                            desc="Identifying jdl and jtl genes: 2 of 2"):
        del_idx_i, del_idx_j = del_idx_pair
        if np.where(jnz_copy == del_idx_j) < np.where(jnz_copy == del_idx_i):
            with model:
                model.genes[del_idx_i].knock_out()
                model.genes[del_idx_j].knock_out()
                sol_ko_ij = model.optimize()
                if sol_ko_ij.objective_value < cutoff * gr_wt and \
                        sol_ko_ij.status != 'infeasible':
                    jdl_genes_idx.append([int(del_idx_i), int(del_idx_j)])

                elif sol_ko_ij.status == 'infeasible':
                    sol_ko_ij = model.optimize()
                    if sol_ko_ij.objective_value < cutoff * gr_wt or \
                            math.isnan(sol_ko_ij.objective_value) is True:
                        jdl_genes_idx.append([int(del_idx_i), int(del_idx_j)])
                        # continue

                    newnnz_ij = np.flatnonzero(sol_ko_ij.fluxes)
                    newnnz_ij_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                                      for rxn_idx in newnnz
                                                      if model.reactions[rxn_idx].gene_reaction_rule != ''}
                    newnnz_ij_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                                   for single_frozenset
                                                                   in newnnz_rxns_genes_frozenset
                                                                   for genes in single_frozenset]))
                    jnz_ij = np.setdiff1d(newnnz_ij_rxns_genes_idx,
                                          jnz_rxns_genes_idx)
                    for del_idx_k in jnz_ij:
                        with model:
                            sol_ko_ijk = model.slim_optimize()
                            if sol_ko_ijk < cutoff * gr_wt or \
                                    math.isnan(sol_ko_ijk) is True:
                                jtl_genes_idx.append([int(del_idx_i),
                                                      int(del_idx_j),
                                                      int(del_idx_k)])

                    for del_idx_k in jnz_copy:
                        with model:
                            if np.where(jnz_copy == del_idx_k) < \
                                    np.where(jnz_copy == del_idx_j):
                                sol_ko_ijk = model.slim_optimize()
                                if sol_ko_ijk < cutoff * gr_wt or \
                                        math.isnan(sol_ko_ijk) is True:
                                    jtl_genes_idx.append([int(del_idx_i),
                                                          int(del_idx_j),
                                                          int(del_idx_k)])

    # Eliminate double lethal reaction deletions in triple lethal reacutions
    jdl_genes_idx = np.array(jdl_genes_idx)
    jtl_genes_idx = np.array(jtl_genes_idx)

    temporary = []
    g = np.zeros(jdl_genes_idx.shape[0])
    for del_idx_i in jtl_genes_idx:
        for del_idx_j in jdl_genes_idx:
            g[np.where(jdl_genes_idx == del_idx_j)[0][0]] = \
                    np.sum(np.in1d(del_idx_i, del_idx_j))
            # if g[np.where(jdl_idx == del_idx_j)[0][0]] >= 2:
            #     break
        if np.max(g) < 2:
            temporary.append(del_idx_i)

    jtl_genes_idx_new = np.array(temporary)
    jtl_genes_idx_final = np.unique(np.sort(jtl_genes_idx_new),
                                    axis=0).tolist()

    # Indices -> Genes
    jdl_genes_idx = jdl_genes_idx.tolist()
    jdl_genes = [model.genes.get_by_any(gene_pair_idx) for gene_pair_idx
                 in jdl_genes_idx]

    jtl_genes = [model.genes.get_by_any(gene_triplet_idx) for gene_triplet_idx
                 in jtl_genes_idx_final]

    return (jsl_genes, jdl_genes, jtl_genes)
