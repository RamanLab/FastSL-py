# -*- coding: utf-8 -*-

import numpy as np
from tqdm import tqdm
import math
from singleSL_genes import singleSL_genes

def doubleSL_genes(model,cutoff,solver):
    '''
    Analysis for double lethal genes
    '''
    model.solver = solver
    solWT = model.optimize()
    grWT = solWT.objective_value

    Jnz_rxns_idx = np.flatnonzero(solWT.fluxes) # gives the reaction indices
    Jnz_rxns_genes_frozenset = {model.reactions[rxn_idx].genes for rxn_idx in Jnz_rxns_idx if model.reactions[rxn_idx].gene_reaction_rule != ''} # gives the set of the frozensets of genes associated with the Jnz reactions; removes duplicates since set
    Jnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes) for single_frozenset in Jnz_rxns_genes_frozenset for genes in single_frozenset])) # gives the indices of the genes obtained; removes duplicates too

    Jsl = singleSL_genes(model=model,cutoff=cutoff,solver=solver)
    Jsl_idx = [model.genes.index(Jsl_id) for Jsl_id in Jsl]

    Jnz_copy = np.setdiff1d(Jnz_rxns_genes_idx,Jsl_idx) # Jnz-Jsl

    Jdl_idx = []

    for delIdx_i in tqdm(Jnz_copy,desc="Identifying Jdl genes - Part 1 of 2"):
        with model:
            model.genes[delIdx_i].knock_out()
            solKO_i = model.optimize()
            newnnz = np.flatnonzero(solKO_i.fluxes)
            newnnz_rxns_genes_frozenset = {model.reactions[rxn_idx].genes for rxn_idx in newnnz if model.reactions[rxn_idx].gene_reaction_rule != ''}
            newnnz_rxns_genes_idx = np.unique(np.array([model.genes.index(genes) for single_frozenset in newnnz_rxns_genes_frozenset for genes in single_frozenset]))
            Jnz_i = np.setdiff1d(newnnz_rxns_genes_idx,Jnz_rxns_genes_idx)

            for delIdx_j in Jnz_i:
                with model:
                    model.genes[delIdx_j].knock_out()
                    solKO_ij = model.slim_optimize()
                    if solKO_ij < cutoff * grWT or math.isnan(solKO_ij) == True:
                        Jdl_idx.append([int(delIdx_i),int(delIdx_j)])

    for delIdx_i in tqdm(Jnz_copy,desc="Identifying Jdl genes - Part 2 of 2"):
        for delIdx_j in Jnz_copy:
            if np.where(Jnz_copy==delIdx_j) < np.where(Jnz_copy==delIdx_i):
                with model:
                    model.genes[delIdx_i].knock_out()
                    model.genes[delIdx_j].knock_out()
                    solKO_ij = model.slim_optimize()
                    if solKO_ij < cutoff * grWT or math.isnan(solKO_ij) == True:
                        Jdl_idx.append([int(delIdx_i),int(delIdx_j)])

    Jdl = [model.genes.get_by_any(gene_pair_idx) for gene_pair_idx in Jdl_idx]
    return (Jsl,Jdl)
