#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Contains functions required for CLI.'''

import logging

import click
from lxml import etree

import cobra
from fastsl.genes import double_genes, single_genes
from fastsl.reactions import double_reactions, single_reactions
from fastsl.util import generate_dir, write_file

LOGGER = logging.getLogger(__name__)


@click.command()
@click.version_option()
@click.option('--cutoff', type=float, default=0.01,
              help='cutoff for the growth rate')
@click.option('--order', type=int, default=2,
              help='order of synthetic lethals')
@click.option('--elilist',
              help='elimination list for model')
@click.option('--atpm', type=str, default='atpm',
              help='ID of ATP maintenance reaction')
@click.option('--solver', type=str, default='glpk_exact',
              help='LP solver')
@click.option('--parallel', is_flag=True, default=False,
              help='parallel version')
@click.option('--genes', is_flag=True, default=False,
              help='fast-sl genes')
@click.option('--gen-elilist', is_flag=True, default=False,
              help='generates the elimination list for the model')
@click.argument('model')
# TODO: refactor branches into functions
def main(model, cutoff, order, elilist, atpm,
         solver, parallel, genes, gen_elilist):
    '''An efficient tool for analysing synthetic lethal reactions/genes
    in genome-wide metabolic networks'''

    logging.basicConfig(filename='fast-sl.log',
                       format='%(asctime)s : %(message)s',
                       datefmt='%m/%d/%Y %I:%M:%S %p',
                       level=logging.INFO)
    # PARSE SBML MODEL
    model = cobra.io.read_sbml_model(model)

    # +---------------------------------------------+
    # | This section generates the elimination list |
    # | for a model.                                |
    # +---------------------------------------------+
    if gen_elilist is True:
        exchange_reactions_objects = model.exchanges
        exchange_reactions_index = [model.reactions.index(reaction_id)
                                    for reaction_id in exchange_reactions_objects]
        exchange_reactions_list = [model.reactions[reaction_idx].id
                                   for reaction_idx in exchange_reactions_index]

        # generate XML
        root = etree.Element('elimination-list')

        for reaction_id in exchange_reactions_list:
            child = etree.SubElement(root, 'reaction-id')
            child.text = reaction_id

        tree = etree.ElementTree(root)

        tree.write('{}_elimination_list.xml'.format(model),
                   encoding='utf-8',
                   xml_declaration=True,
                   pretty_print=True)

        LOGGER.info('%s_elimination_list.xml generated.', model)

    # +------------------------------------------+
    # | This section performs FastSL obeying the |
    # | configuration provided.                  |
    # +------------------------------------------+
    else:
        # CREATE DIRECTORY BASED ON WHETHER INPUT IS REACTIONS/GENES
        if genes is True:
            results_dir_path = generate_dir(model.id, "genes")
        else:
            results_dir_path = generate_dir(model.id, "reactions")

            # ELIMINATION LIST DATA PROCESSING
            if not elilist:
                # error handling for adding ATP maintenance reaction
                # in the elimination list if no elilist provided
                try:
                    elilist_data = [model.reactions.get_by_id(atpm).id]
                    LOGGER.info(('%s found in model and added to elimination '
                                 'list.'), atpm)
                except KeyError:
                    LOGGER.critical(('%s not found in model. Provide a valid '
                                     'ATP maintenance reaction ID or provide '
                                     'an elimination list.'), atpm)
                    raise SystemExit('%s not found in model. Provide a valid '
                                     'ATP maintenance reaction ID or provide '
                                     'an elimination list.' % atpm)

            else:
                # elimination list parsing
                elilist_tree = etree.parse(elilist)
                elilist_data = [data.text for data in
                                elilist_tree.iter(tag='reaction-id')]

        # ======= SERIAL/PARALLEL HANDLING FOR REACTIONS =======
        if parallel is False and genes is False:

            if order == 1:
                jsl = single_reactions(model, elilist_data, cutoff,
                                       solver=solver, cpu_count=1)
                write_file(results_dir_path, model.id,
                           "reactions", "single", jsl)

                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info('Single lethal reactions: %s', len(jsl))


            elif order == 2:
                jsl, jdl = double_reactions(model, elilist_data, cutoff,
                                            solver=solver)

                write_file(results_dir_path, model.id,
                           "reactions", "single", jsl)
                write_file(results_dir_path, model.id,
                           "reactions", "double", jdl)

                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info(('Single lethal reactions: %s\t'
                             'Double lethal reactions: %s'),
                            len(jsl), len(jdl))

        elif parallel is True and genes is False:

            if order == 1:
                jsl = single_reactions(model, elilist_data, cutoff,
                                       solver=solver)
                write_file(results_dir_path,
                           model.id,
                           "reactions",
                           "single",
                           jsl)
                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info('Single lethal reactions: %s', len(jsl))

            elif order == 2:
                jsl, jdl = double_reactions(model, elilist_data, cutoff,
                                            solver=solver)

                write_file(results_dir_path, model.id,
                           "reactions", "single", jsl)
                write_file(results_dir_path, model.id,
                           "reactions", "double", jdl)

                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info(('Single lethal reactions: %s\t'
                             'Double lethal reactions: %s'),
                            len(jsl), len(jdl))

        # ======= SERIAL/PARALLEL HANDLING FOR GENES =======
        if parallel is False and genes is True:

            if order == 1:
                jsl = single_genes(model, cutoff, solver=solver)

                write_file(results_dir_path, model.id,
                           "genes", "single", jsl)

                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info('Single lethal genes: %s', len(jsl))

            elif order == 2:
                jsl, jdl = double_genes(model, cutoff, solver=solver)

                write_file(results_dir_path, model.id,
                           "genes", "single", jsl)
                write_file(results_dir_path, model.id,
                           "genes", "double", jdl)

                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info(('Single lethal genes: %s\t'
                             'Double lethal genes: %s'),
                            len(jsl), len(jdl))

        elif parallel is True and genes is True:

            if order == 1:
                jsl = single_genes(model, cutoff, solver=solver)

                write_file(results_dir_path, model.id,
                           "genes", "single", jsl)

                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info('Single lethal genes: %s', len(jsl))

            elif order == 2:
                jsl, jdl = double_genes(model, cutoff, solver=solver)

                write_file(results_dir_path, model.id,
                           "genes", "single", jsl)
                write_file(results_dir_path, model.id,
                           "genes", "double", jdl)

                LOGGER.info('========== %s ==========', model.id)
                LOGGER.info(('Single lethal genes: %s\t'
                             'Double lethal genes: %s'),
                            len(jsl), len(jdl))


if __name__ == '__main__':
    main()
