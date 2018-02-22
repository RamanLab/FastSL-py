#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import argparse, logging, cobra
from lxml import etree
from singleSL_genes import singleSL_genes
from doubleSL_genes import doubleSL_genes
# from tripleSL import tripleSL

def fastSL_genes(model,cutoff,order,solver):
    '''Identifies synthetic lethal genes based on the parameters provided'''

    # model parsing
    model = cobra.io.read_sbml_model(model)

   # order handling
    if order == 1:
        Jsl = singleSL_genes(model,cutoff,solver)
        logging.info('Jsl: {}'.format(len(Jsl)))
    elif order == 2:
        Jsl,Jdl = doubleSL_genes(model,cutoff,solver)
        logging.info('Jsl: {}; Jdl: {}'.format(len(Jsl), len(Jdl)))
    # elif order == 3:
    #     Jsl,Jdl,Jtl = tripleSL(model,cutoff,eliListData,atpm,solver)
        # logging.info('Jsl: {}; Jdl: {}; Jtl: {}'.format(len(Jsl), len(Jdl), len(Jtl)))

def main():
    # command-line argument parser
    parser = argparse.ArgumentParser(description='Identification of synthetic lethal genes of a genome-wide metabolic network model')
    parser.add_argument('model',nargs='?',type=str,help='model .xml file')
    parser.add_argument('--cutoff',nargs='?',default=0.01,type=int,help='cut-off value')
    parser.add_argument('--order',nargs='?',default=2,type=int,help='order')
    parser.add_argument('--solver',nargs='?',default='glpk',type=str,help='solver used')

    args = parser.parse_args()

    # Logger initiation
    logging.getLogger().setLevel(logging.INFO)

    # FastSL
    fastSL_genes(model=args.model,cutoff=args.cutoff,order=args.order,solver=args.solver)

if __name__ == '__main__':
    main()
