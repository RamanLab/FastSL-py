# -*- coding: utf-8 -*-

import argparse
import logging

import cobra
from lxml import etree

from fastSL_genes import singleSL_genes, doubleSL_genes, tripleSL_genes


def fastSL_genes(model, cutoff, order, solver):
    '''Identifies synthetic lethal reactions'''

    # model parsing
    model = cobra.io.read_sbml_model(model)

    # order handling
    if order == 1:
        Jsl = singleSL_genes(model,
                             cutoff,
                             solver)
        logging.info('Jsl:{}\n{}'
                     .format(len(Jsl),
                             Jsl))
    elif order == 2:
        Jsl, Jdl = doubleSL_genes(model,
                                  cutoff,
                                  solver)
        logging.info('Jsl:{}\n{};\nJdl:{}\n{}'
                     .format(len(Jsl),
                             Jsl,
                             len(Jdl),
                             Jdl))
    elif order == 3:
        Jsl, Jdl, Jtl = tripleSL_genes(model,
                                       cutoff,
                                       solver)
        logging.info('Jsl:{}\n{};\nJdl:{}\n{};\nJtl:{}\n{}'
                     .format(len(Jsl),
                             Jsl,
                             len(Jdl),
                             Jdl,
                             len(Jtl),
                             len(Jtl)))


def main():
    # command-line argument parser
    parser = argparse.ArgumentParser(description='Identification of synthetic\
                                                  lethal genes of a\
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
    parser.add_argument('--solver',
                        nargs='?',
                        default='glpk',
                        type=str,
                        help='solver used')

    args = parser.parse_args()

    # Logger initiation
    logging.getLogger().setLevel(logging.INFO)

    # FastSL
    fastSL_genes(args.model,
                 args.cutoff,
                 args.order,
                 args.solver)

if __name__ == '__main__':
    main()