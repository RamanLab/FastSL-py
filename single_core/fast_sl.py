# -*- coding: utf-8 -*-

import argparse
import logging

import cobra
from lxml import etree

def fast_sl(model, cutoff, order, atpm, solver, elilist=None):
    """Identify synthetic lethal reactions."""
    # from lxml import etree
    from rxns import single_sl, double_sl, triple_sl

    # model parsing
    model = cobra.io.read_sbml_model(model)

    # elilist data processing
    if not elilist:
        # error handling for adding ATP maintenance reaction
        # in the elimination list if no elilist provided
        try:
            elilist_data = model.reactions.get_by_id(atpm).id
            logging.info('%s found in model and added to elimination list.',
                         atpm)
        except KeyError:
            logging.info('%s not found in the model.', atpm)
    elif elilist:
        # elimination list parsing
        elilist_tree = etree.parse(elilist)
        elilist_data = [data.text for data in
                        elilist_tree.iter(tag='reaction-id')]

    # order handling
    if order == 1:
        jsl = single_sl(model,
                        cutoff,
                        elilist_data,
                        solver)
        logging.info('jsl:%s\n%s', len(jsl), jsl)
    elif order == 2:
        jsl, jdl = double_sl(model,
                             cutoff,
                             elilist_data,
                             solver)
        logging.info('jsl:%s\n%s;\njdl:%s\n%s',
                     len(jsl), jsl,
                     len(jdl), jdl)
    elif order == 3:
        jsl, jdl, jtl = triple_sl(model,
                                  cutoff,
                                  elilist_data,
                                  solver)
        logging.info('jsl:%s\n%s;\njdl:%s\n%s;\njtl:%s\n%s',
                     len(jsl), jsl,
                     len(jdl), jdl,
                     len(jtl), jtl)


def fast_sl_genes(model, cutoff, order, solver):
    """Identify synthetic lethal genes."""
    from genes import single_sl_genes, double_sl_genes, triple_sl_genes

    # model parsing
    model = cobra.io.read_sbml_model(model)

    # order handling
    if order == 1:
        jsl = single_sl_genes(model,
                              cutoff,
                              solver)
        logging.info('jsl:%s\n%s', len(jsl), jsl)
    elif order == 2:
        jsl, jdl = double_sl_genes(model,
                                   cutoff,
                                   solver)
        logging.info('jsl:%s\n%s;\njdl:%s\n%s',
                     len(jsl), jsl,
                     len(jdl), jdl)
    elif order == 3:
        jsl, jdl, jtl = triple_sl_genes(model,
                                        cutoff,
                                        solver)
        logging.info('jsl:%s\n%s;\njdl:%s\n%s;\njtl:%s\n%s',
                     len(jsl), jsl,
                     len(jdl), jdl,
                     len(jtl), jtl)


def main():
    """Take command line arguments and run the necessary functions."""
    # command-line argument parser
    parser = argparse.ArgumentParser(description='Identification of synthetic\
                                                  lethal reactions/genes of a\
                                                  genome-wide metabolic\
                                                  network model')
    parser.add_argument('model',
                        nargs='?',
                        type=str,
                        help='model .xml file')
    parser.add_argument('--cutoff',
                        nargs='?',
                        default=0.01,
                        type=int,
                        help='cut-off value')
    parser.add_argument('--order',
                        nargs='?',
                        default=2,
                        type=int,
                        help='order')
    parser.add_argument('--elilist',
                        nargs='?',
                        type=str,
                        help='.xml file of elimination reactions')
    parser.add_argument('--atpm',
                        nargs='?',
                        default='ATPM',
                        type=str,
                        help='ID of ATP maintenance reaction')
    parser.add_argument('--solver',
                        nargs='?',
                        default='glpk_exact',
                        type=str,
                        help='solver used')
    parser.add_argument('--genes',
                        nargs='?',
                        default=0,
                        type=int,
                        help='0: reactions, 1: genes')
    args = parser.parse_args()

    # Logger initiation
    logging.getLogger().setLevel(logging.INFO)

    if args.genes == 0:
        # Fast_sl
        fast_sl(args.model,
                args.cutoff,
                args.order,
                args.atpm,
                args.solver,
                args.elilist)
    else:
        # Fast_sl_genes
        fast_sl_genes(args.model,
                      args.cutoff,
                      args.order,
                      args.solver)


if __name__ == '__main__':
    main()
