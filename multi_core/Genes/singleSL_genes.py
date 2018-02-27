# -*- coding: utf-8 -*-

import numpy as np
from tqdm import tqdm
import math
from multiprocessing import Pool, cpu_count


def _filter_single_lethal_genes(model, cutoff, grWT, delIdx_i):
    with model:
        model.genes[delIdx_i].knock_out()
        solKO_i = model.slim_optimize()
        if solKO_i < cutoff * grWT or math.isnan(solKO_i) is True:
            return int(delIdx_i)


def _filter_single_lethal_genes_worker(genes_idx):
    global _model, _cutoff, _grWT
    return _filter_single_lethal_genes(_model, _cutoff, _grWT, genes_idx)


def _init_worker(model, cutoff, grWT):
    global _model, _cutoff, _grWT
    _model = model
    _cutoff = cutoff
    _grWT = grWT


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
    # Identify Single Lethal Gene Deletions...

    # Multiprocessing
    chunk_size = Jnz_rxns_genes_idx.shape[0] // cpu_count()
    processes = cpu_count()
    # initiate pool workers
    pool = Pool(processes,
                initializer=_init_worker,
                initargs=(model, cutoff, grWT))
    Jsl_genes_idx = list(
                         filter(
                                lambda rxn_idx: rxn_idx is not None,
                                pool.imap_unordered(
                                                    _filter_single_lethal_genes_worker,
                                                    Jnz_rxns_genes_idx,
                                                    chunksize=chunk_size)))
    pool.close()  # stop pool workers
    pool.join()  # terminate pool workers

    # for delIdx_i in tqdm(Jnz_rxns_genes_idx,desc="Identifying Jsl genes"):
    #     with model:
    #         model.genes[delIdx_i].knock_out()
    #         solKO_i = model.slim_optimize()
    #         if solKO_i < cutoff * grWT or math.isnan(solKO_i) == True:
    #             Jsl_genes_idx.append(int(delIdx_i))

    Jsl = model.genes.get_by_any(Jsl_genes_idx)
    return Jsl
