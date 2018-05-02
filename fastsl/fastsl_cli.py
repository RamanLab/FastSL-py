#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import sys
from pathlib import Path
import csv

import click
import cobra

from lxml import etree

LOGGER = logging.getLogger(__name__)


def generate_dir(model_name, of_type):
    '''Generate the required directory for output.'''
    dir_path = Path.cwd() / 'fast-sl-results/{}/{}'.format(model_name,
                                                           of_type)
    dir_path.mkdir(parents=True, exist_ok=True)

    logging.info('%s generated.', str(dir_path.resolve()))
    return dir_path


def write_file(in_dir, model_name, of_type, lethality_stage, data):
    '''Writes the output file for the operation.'''
    file_path = in_dir / '{}_{}_lethal_{}.csv'.format(model_name,
                                                      lethality_stage,
                                                      of_type)

    file_path.touch()
    logging.info('%s generated.', str(file_path.resolve()))

    with file_path.open('w') as file:
        if lethality_stage == 'single':
            csv.writer(file).writerow(data)
        elif lethality_stage == 'double' or lethality_stage == 'triple':
            csv.writer(file).writerows(data)


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
def main(model, cutoff, order, elilist, atpm, solver, parallel, genes,
         gen_elilist):
    '''An efficient tool for analysing synthetic lethal reactions/genes
    in genome-wide metabolic networks'''
    logging.basicConfig(filename='fast-sl.log',
                        format='%(asctime)s : %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)
    # model parsing
    model = cobra.io.read_sbml_model(model)

    # performs elilist generation
    if gen_elilist is True:
        exchange_reactions_objects = model.exchanges
        exchange_reactions_index = [model.reactions.index(reaction_id)
                                    for reaction_id in exchange_reactions_objects]
        exchange_reactions_list = [model.reactions[reaction_idx].id
                                   for reaction_idx in exchange_reactions_index]

        # xml generation
        root = etree.Element('elimination-list')

        for reaction_id in exchange_reactions_list:
            child = etree.SubElement(root, 'reaction-id')
            child.text = reaction_id

        tree = etree.ElementTree(root)

        tree.write('{}_elimination_list.xml'.format(model),
                   encoding='utf-8',
                   xml_declaration=True,
                   pretty_print=True)

        logging.info('%s_elimination_list.xml generated.', model)
    # performs fast-sl
    else:
        # directory creation based on whether input is reactions/genes
        if genes is True:
            results_dir_path = generate_dir(model.id, "genes")
        else:
            results_dir_path = generate_dir(model.id, "reactions")
            # elilist data processing
            if not elilist:
                # error handling for adding ATP maintenance reaction
                # in the elimination list if no elilist provided
                try:
                    elilist_data = [model.reactions.get_by_id(atpm).id]
                    logging.info(('%s found in model and added to elimination '
                                  'list.'),
                                 atpm)
                except KeyError:
                    logging.critical(('%s not found in model. Enter a valid ATP '
                                      'maintenance reaction ID or provide an '
                                      'elimination list.'), atpm)
                    sys.exit()
            elif elilist:
                # elimination list parsing
                elilist_tree = etree.parse(elilist)
                elilist_data = [data.text for data in
                                elilist_tree.iter(tag='reaction-id')]

        # serial/parallel handling for reactions
        if parallel is False and genes is False:
            from fastsl.rxns import single_sl, double_sl
            # Order 1
            if order == 1:
                jsl = single_sl(model,
                                cutoff,
                                elilist_data,
                                solver)
                write_file(results_dir_path,
                           model.id,
                           "reactions",
                           "single",
                           jsl)
                logging.info('========== %s ==========', model.id)
                logging.info('Jsl:%s', len(jsl))
            # Order 2
            elif order == 2:
                jsl, jdl = double_sl(model,
                                     cutoff,
                                     elilist_data,
                                     solver)
                write_file(results_dir_path,
                           model.id,
                           "reactions",
                           "single",
                           jsl)
                write_file(results_dir_path,
                           model.id,
                           "reactions",
                           "double",
                           jdl)
                logging.info('========== %s ==========', model.id)
                logging.info('Jsl:%s\tJdl:%s', len(jsl), len(jdl))

        elif parallel is True and genes is False:
            from fastsl.parallel_rxns import (
                parallel_single_sl,
                parallel_double_sl)
            # Parallel order 1
            if order == 1:
                jsl = parallel_single_sl(model,
                                         cutoff,
                                         elilist_data,
                                         solver)
                write_file(results_dir_path,
                           model.id,
                           "reactions",
                           "single",
                           jsl)
                logging.info('========== %s ==========', model.id)
                logging.info('Jsl:%s', len(jsl))
            # Parallel order 2
            elif order == 2:
                jsl, jdl = parallel_double_sl(model,
                                              cutoff,
                                              elilist_data,
                                              solver)
                write_file(results_dir_path,
                           model.id,
                           "reactions",
                           "single",
                           jsl)
                write_file(results_dir_path,
                           model.id,
                           "reactions",
                           "double",
                           jdl)
                logging.info('========== %s ==========', model.id)
                logging.info('Jsl:%s\tJdl:%s',
                             len(jsl),
                             len(jdl))

        # serial/parallel handling for genes
        if parallel is False and genes is True:
            from fastsl.genes import single_sl_genes, double_sl_genes
            # Order 1
            if order == 1:
                jsl = single_sl_genes(model,
                                      cutoff,
                                      solver)
                write_file(results_dir_path,
                           model.id,
                           "genes",
                           "single",
                           jsl)
                logging.info('========== %s ==========', model.id)
                logging.info('jsl genes:%s', len(jsl))
            # Order 2
            elif order == 2:
                jsl, jdl = double_sl_genes(model,
                                           cutoff,
                                           solver)
                write_file(results_dir_path,
                           model.id,
                           "genes",
                           "single",
                           jsl)
                write_file(results_dir_path,
                           model.id,
                           "genes",
                           "double",
                           jdl)
                logging.info('========== %s ==========', model.id)
                logging.info('Jsl genes:%s\tJdl genes:%s', len(jsl), len(jdl))

        elif parallel is True and genes is True:
            from fastsl.parallel_genes import (
                parallel_single_sl_genes,
                parallel_double_sl_genes)
            # Parallel order 1
            if order == 1:
                jsl = parallel_single_sl_genes(model,
                                               cutoff,
                                               solver)
                write_file(results_dir_path,
                           model.id,
                           "genes",
                           "single",
                           jsl)
                logging.info('========== %s ==========', model.id)
                logging.info('Jsl genes:%s', len(jsl))
            # Parallel order 2
            elif order == 2:
                jsl, jdl = parallel_double_sl_genes(model,
                                                    cutoff,
                                                    solver)
                write_file(results_dir_path,
                           model.id,
                           "genes",
                           "single",
                           jsl)
                write_file(results_dir_path,
                           model.id,
                           "genes",
                           "double",
                           jdl)
                logging.info('========== %s ==========', model.id)
                logging.info('Jsl genes:%s\tJdl genes:%s',
                             len(jsl),
                             len(jdl))


if __name__ == '__main__':
    main()
