# -*- coding: utf-8 -*-

import math
from os import cpu_count

import numpy as np
from joblib import Parallel, delayed


def _single_lethal_reactions(model, cutoff, grWT, delIdx_i):
    with model:
        model.reactions[delIdx_i].knock_out()
        solKO_i = model.slim_optimize()
        if solKO_i < cutoff * grWT or math.isnan(solKO_i) is True:
            return int(delIdx_i)
        else:
            return None


def singleSL(model, cutoff, eliList, solver):
    '''
    Analysis for single lethal reactions
    '''

    model.solver = solver  # Basic solver configuration
    solWT = model.optimize()  # Identify minNorm flux distribution
    grWT = solWT.objective_value

    # Indices of non-zero flux reactions including exchange/diffusion reactions
    Jnz_before_filtering = np.flatnonzero(solWT.fluxes)

    # Indices of exchange/diffusion reactions (not considered)
    eliIdx = [model.reactions.index(reaction_id) for reaction_id in eliList]

    Jnz = np.setdiff1d(Jnz_before_filtering, eliIdx)  # Jnz

    # Identify Single Lethal Reaction Deletions

    chunk_size = Jnz.shape[0] // cpu_count()  # integer division

    Jsl_idx = list(
                   filter(
                          lambda rxn_idx: rxn_idx is not None,
                          Parallel(n_jobs=cpu_count(),
                                   # threading performs better than
                                   # multiprocessing in only deletions
                                   backend='threading',
                                   verbose=3,
                                   batch_size=chunk_size)(
                                   delayed(_single_lethal_reactions)
                                   (model,
                                    cutoff,
                                    grWT,
                                    delIdx_i) for delIdx_i in Jnz)))

    # Indices -> Reaction objects
    Jsl = model.reactions.get_by_any(Jsl_idx)

    return Jsl
