# -*- coding: utf-8 -*-

import argparse
import logging

import cobra
from lxml import etree

from fastSL_rxns import singleSL, doubleSL, tripleSL


def fastSL(model, cutoff, order, atpm, solver, eliList=None):
    '''
    Identifies synthetic lethal reactions
    '''

    # model parsing
    model = cobra.io.read_sbml_model(model)

    # eliList data processing
    if not eliList:
        # error handling for adding ATP maintenance reaction
        # in the elimination list if no eliList provided
        try:
            eliListData = model.reactions.get_by_id(atpm).id
            logging.info('%s found in model and added to elimination list.',
                         atpm)
        except KeyError:
            logging.info('%s not found in the model.', atpm)
    elif eliList:
        # elimination list parsing
        eliListTree = etree.parse(eliList)
        eliListData = [data.text for data in
                       eliListTree.iter(tag='reaction-id')]

    # order handling
    if order == 1:
        Jsl = singleSL(model,
                       cutoff,
                       eliListData,
                       solver)
        logging.info('Jsl:%s\n%s', len(Jsl), Jsl)
    elif order == 2:
        Jsl, Jdl = doubleSL(model,
                            cutoff,
                            eliListData,
                            solver)
        logging.info('Jsl:%s\n%s;\nJdl:%s\n%s',
                     len(Jsl), Jsl,
                     len(Jdl), Jdl)
    elif order == 3:
        Jsl, Jdl, Jtl = tripleSL(model,
                                 cutoff,
                                 eliListData,
                                 solver)
        logging.info('Jsl:%s\n%s;\nJdl:%s\n%s;\nJtl:%s\n%s',
                     len(Jsl), Jsl,
                     len(Jdl), Jdl,
                     len(Jtl), Jtl)


def main():
    '''
    Takes command line arguments and runs the necessary functions
    '''
    # command-line argument parser
    parser = argparse.ArgumentParser(description='Identification of synthetic\
                                                  lethal reactions of a\
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
    parser.add_argument('--eliList',
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

    args = parser.parse_args()

    # Logger initiation
    logging.getLogger().setLevel(logging.INFO)

    # FastSL
    fastSL(args.model,
           args.cutoff,
           args.order,
           args.atpm,
           args.solver,
           args.eliList)


if __name__ == '__main__':
    main()
