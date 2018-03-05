# -*- coding: utf-8 -*-

import math
from itertools import product

import numpy as np
from tqdm import tqdm


def singleSL_genes(model, cutoff, solver):
    '''
    Analysis for single lethal genes
    '''

    model.solver = solver  # Basic solver configuration
    solWT = model.optimize()  # Identify minNorm flux distribution
    grWT = solWT.objective_value

    # gives the non-zero flux reaction indices
    Jnz_rxns = np.flatnonzero(solWT.fluxes)
    # set of the frozensets of genes associated with the Jnz reactions
    # removes duplicates since set
    Jnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                          for rxn_idx in Jnz_rxns
                          if model.reactions[rxn_idx].gene_reaction_rule != ''}
    # indices of the genes obtained; removes duplicates too
    Jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             Jnz_rxns_genes_set
                                             for genes in single_frozenset]))

    # Identify Single Lethal Gene Deletions

    Jsl_genes_idx = []

    for delIdx_i in tqdm(Jnz_rxns_genes_idx, desc="Identifying Jsl genes"):
        with model:
            model.genes[delIdx_i].knock_out()
            solKO_i = model.slim_optimize()
            if solKO_i < cutoff * grWT or math.isnan(solKO_i) is True:
                Jsl_genes_idx.append(int(delIdx_i))

    # Indices -> Genes
    Jsl_genes = model.genes.get_by_any(Jsl_genes_idx)

    return Jsl_genes


def doubleSL_genes(model, cutoff, solver):
    '''
    Analysis for double lethal genes
    '''

    model.solver = solver  # Basic solver configuration
    solWT = model.optimize()  # Identify minNorm flux distribution
    grWT = solWT.objective_value

    # gives the non-zero flux reaction indices
    Jnz_rxns = np.flatnonzero(solWT.fluxes)
    # set of the frozensets of genes associated with the Jnz reactions
    # removes duplicates and non-gene associated reactions
    Jnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                          for rxn_idx in Jnz_rxns
                          if model.reactions[rxn_idx].gene_reaction_rule != ''}
    # indices of the genes obtained; removes duplicates too
    Jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             Jnz_rxns_genes_set
                                             for genes in single_frozenset]))

    # Identify single lethal genes
    Jsl_genes = singleSL_genes(model,
                               cutoff,
                               solver)

    # Indices of single lethal genes
    Jsl_genes_idx = [model.genes.index(Jsl_id) for Jsl_id in Jsl_genes]

    Jnz_copy = np.setdiff1d(Jnz_rxns_genes_idx, Jsl_genes_idx)  # Jnz-Jsl

    # Makes rxn pairs to remove nested for-loops in phase 2
    Jnz_copy_phase_2 = [gene_pair for gene_pair in product(Jnz_copy, repeat=2)]

    # Identify double lethal reactions

    Jdl_genes_idx = []

    # Phase 1:
    for delIdx_i in tqdm(Jnz_copy, desc="Identifying Jdl genes: 1 of 2"):
        with model:
            model.genes[delIdx_i].knock_out()
            solKO_i = model.optimize()
            newnnz = np.flatnonzero(solKO_i.fluxes)
            newnnz_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                           for rxn_idx in newnnz
                                           if model.reactions[rxn_idx].gene_reaction_rule != ''}
            newnnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                        for single_frozenset
                                                        in newnnz_rxns_genes_frozenset
                                                        for genes in single_frozenset]))
            Jnz_i = np.setdiff1d(newnnz_rxns_genes_idx,
                                 Jnz_rxns_genes_idx)

            for delIdx_j in Jnz_i:
                with model:
                    model.genes[delIdx_j].knock_out()
                    solKO_ij = model.slim_optimize()
                    if solKO_ij < cutoff * grWT or \
                            math.isnan(solKO_ij) is True:
                        Jdl_genes_idx.append([int(delIdx_i),
                                              int(delIdx_j)])

    # Phase 2:
    for delIdx_pair in tqdm(Jnz_copy_phase_2,
                            desc="Identifying Jdl genes: 2 of 2"):
        delIdx_i, delIdx_j = delIdx_pair
        if np.where(Jnz_copy == delIdx_j) < np.where(Jnz_copy == delIdx_i):
            with model:
                model.genes[delIdx_i].knock_out()
                model.genes[delIdx_j].knock_out()
                solKO_ij = model.slim_optimize()
                if solKO_ij < cutoff * grWT or \
                        math.isnan(solKO_ij) is True:
                    Jdl_genes_idx.append([int(delIdx_i), int(delIdx_j)])

    # Indices -> Genes
    Jdl_genes = [model.genes.get_by_any(gene_pair_idx) for gene_pair_idx
                 in Jdl_genes_idx]

    return (Jsl_genes, Jdl_genes)


def tripleSL_genes(model, cutoff, solver):
    '''
    Analysis for triple lethal reactions
    '''

    model.solver = solver  # Basic solver configuration
    solWT = model.optimize()  # Identify minNorm flux distribution
    grWT = solWT.objective_value

     # gives the non-zero flux reaction indices
    Jnz_rxns = np.flatnonzero(solWT.fluxes)
    # set of the frozensets of genes associated with the Jnz reactions
    # removes duplicates and non-gene associated reactions
    Jnz_rxns_genes_set = {model.reactions[rxn_idx].genes
                          for rxn_idx in Jnz_rxns
                          if model.reactions[rxn_idx].gene_reaction_rule != ''}
    # indices of the genes obtained; removes duplicates too
    Jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                             for single_frozenset in
                                             Jnz_rxns_genes_set
                                             for genes in single_frozenset]))

    # Identify single lethal genes
    Jsl_genes = singleSL_genes(model,
                               cutoff,
                               solver)

    # Indices of single lethal genes
    Jsl_genes_idx = [model.genes.index(Jsl_id) for Jsl_id in Jsl_genes]

    Jnz_copy = np.setdiff1d(Jnz_rxns_genes_idx, Jsl_genes_idx)  # Jnz-Jsl

    # Makes rxn pairs to remove nested for-loops in phase 2
    Jnz_copy_phase_2 = [gene_pair for gene_pair in product(Jnz_copy, repeat=2)]

    # Identify triple lethal reactions

    Jdl_genes_idx = []
    Jtl_genes_idx = []

    # Phase 1:
    for delIdx_i in tqdm(Jnz_copy,
                         desc="Identifying Jdl and Jtl genes: 1 of 2"):
        with model:
            model.genes[delIdx_i].knock_out()
            solKO_i = model.optimize()
            newnnz = np.flatnonzero(solKO_i.fluxes)
            newnnz_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                           for rxn_idx in newnnz
                                           if model.reactions[rxn_idx].gene_reaction_rule != ''}
            newnnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                        for single_frozenset
                                                        in newnnz_rxns_genes_frozenset
                                                        for genes in single_frozenset]))
            Jnz_i = np.setdiff1d(newnnz_rxns_genes_idx,
                                 Jnz_rxns_genes_idx)

            for delIdx_j in Jnz_i:
                with model:
                    model.genes[delIdx_j].knock_out()
                    solKO_ij = model.optimize()
                    if solKO_ij.objective_value < cutoff * grWT or \
                            solKO_ij.status != 'infeasible':
                        Jdl_genes_idx.append([int(delIdx_i),
                                              int(delIdx_j)])

                    elif solKO_ij.status == 'infeasible':
                        solKO_ij = model.optimize()
                        if solKO_ij.objective_value < cutoff * grWT and \
                                math.isnan(solKO_ij.objective_value) is True:
                            Jdl_genes_idx.append([int(delIdx_i),
                                                  int(delIdx_j)])
                            # continue

                        newnnz_ij = np.flatnonzero(solKO_ij.fluxes)
                        newnnz_ij_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                                          for rxn_idx in newnnz
                                                          if model.reactions[rxn_idx].gene_reaction_rule != ''}
                        newnnz_ij_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                                       for single_frozenset
                                                                       in newnnz_rxns_genes_frozenset
                                                                       for genes in single_frozenset]))
                        Jnz_ij = np.setdiff1d(newnnz_ij_rxns_genes_idx,
                                              Jnz_rxns_genes_idx)
                        for delIdx_k in Jnz_ij:
                            with model:
                                model.genes[delIdx_k].knock_out()
                                solKO_ijk = model.slim_optimize()
                                if solKO_ijk < cutoff * grWT or \
                                        math.isnan(solKO_ijk) is True:
                                    Jtl_genes_idx.append([int(delIdx_i),
                                                          int(delIdx_j),
                                                          int(delIdx_k)])

    # Phase 2:
    for delIdx_pair in tqdm(Jnz_copy_phase_2,
                            desc="Identifying Jdl and Jtl genes: 2 of 2"):
        delIdx_i, delIdx_j = delIdx_pair
        if np.where(Jnz_copy == delIdx_j) < np.where(Jnz_copy == delIdx_i):
            with model:
                model.genes[delIdx_i].knock_out()
                model.genes[delIdx_j].knock_out()
                solKO_ij = model.optimize()
                if solKO_ij.objective_value < cutoff * grWT and \
                        solKO_ij.status != 'infeasible':
                    Jdl_genes_idx.append([int(delIdx_i), int(delIdx_j)])

                elif solKO_ij.status == 'infeasible':
                    solKO_ij = model.optimize()
                    if solKO_ij.objective_value < cutoff * grWT or \
                            math.isnan(solKO_ij.objective_value) is True:
                        Jdl_genes_idx.append([int(delIdx_i), int(delIdx_j)])
                        # continue

                    newnnz_ij = np.flatnonzero(solKO_ij.fluxes)
                    newnnz_ij_rxns_genes_frozenset = {model.reactions[rxn_idx].genes
                                                      for rxn_idx in newnnz
                                                      if model.reactions[rxn_idx].gene_reaction_rule != ''}
                    newnnz_ij_rxns_genes_idx = np.unique(np.array([model.genes.index(genes)
                                                                   for single_frozenset
                                                                   in newnnz_rxns_genes_frozenset
                                                                   for genes in single_frozenset]))
                    Jnz_ij = np.setdiff1d(newnnz_ij_rxns_genes_idx,
                                          Jnz_rxns_genes_idx)
                    for delIdx_k in Jnz_ij:
                        with model:
                            solKO_ijk = model.slim_optimize()
                            if solKO_ijk < cutoff * grWT or \
                                    math.isnan(solKO_ijk) is True:
                                Jtl_genes_idx.append([int(delIdx_i),
                                                int(delIdx_j),
                                                int(delIdx_k)])

                    for delIdx_k in Jnz_copy:
                        with model:
                            if np.where(Jnz_copy == delIdx_k) < \
                                    np.where(Jnz_copy == delIdx_j):
                                solKO_ijk = model.slim_optimize()
                                if solKO_ijk < cutoff * grWT or \
                                        math.isnan(solKO_ijk) is True:
                                    Jtl_genes_idx.append([int(delIdx_i),
                                                          int(delIdx_j),
                                                          int(delIdx_k)])

    # Eliminate double lethal reaction deletions in triple lethal reacutions
    Jdl_genes_idx = np.array(Jdl_genes_idx)
    Jtl_genes_idx = np.array(Jtl_genes_idx)

    temporary = []
    g = np.zeros(Jdl_genes_idx.shape[0])
    for delIdx_i in Jtl_genes_idx:
        for delIdx_j in Jdl_genes_idx:
            g[np.where(Jdl_genes_idx == delIdx_j)[0][0]] = \
                    np.sum(np.in1d(delIdx_i, delIdx_j))
            # if g[np.where(Jdl_idx == delIdx_j)[0][0]] >= 2:
            #     break
        if np.max(g) < 2:
            temporary.append(delIdx_i)

    Jtl_genes_idx_new = np.array(temporary)
    Jtl_genes_idx_final = np.unique(np.sort(Jtl_genes_idx_new), axis=0).tolist()

    # Indices -> Genes
    Jdl_genes_idx = Jdl_genes_idx.tolist()
    Jdl_genes = [model.genes.get_by_any(gene_pair_idx) for gene_pair_idx
                 in Jdl_genes_idx]

    Jtl_genes = [model.genes.get_by_any(gene_triplet_idx) for gene_triplet_idx
                 in Jtl_genes_idx_final]

    return (Jsl_genes, Jdl_genes, Jtl_genes)
