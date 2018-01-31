#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import argparse, logging, cobra
from lxml import etree
from singleSL import singleSL
from doubleSL import doubleSL
from tripleSL import tripleSL

def fastSL(model,cutoff,order,atpm,solver,eliList=None):
    '''Identifies synthetic lethal reactions based on the parameters provided'''

    # model parsing
    model = cobra.io.read_sbml_model(model)

    if not eliList:
        # error handling for adding ATP maintenance reaction in the elimination list
        try:
            eliListData = []
            eliListData.append(model.reactions.get_by_id(atpm).id)
            logging.info('{} found in model and added to elimination list.'.format(atpm))
        except:
            logging.info('{} not found in the model.'.format(atpm))
    elif eliList:
        # elimination list parsing
        eliListTree = etree.parse(eliList)
        eliListData = [data.text for data in eliListTree.iter(tag='reaction-id')]

   # order handling
    if order == 1:
        Jsl = singleSL(model,cutoff,eliListData,atpm,solver)
        logging.info('Jsl: {}'.format(len(Jsl)))
    elif order == 2:
        Jsl,Jdl = doubleSL(model,cutoff,eliListData,atpm,solver)
        logging.info('Jsl: {}; Jdl: {}'.format(len(Jsl), len(Jdl)))
    elif order == 3:
        Jsl,Jdl,Jtl = tripleSL(model,cutoff,eliListData,atpm,solver)
        logging.info('Jsl: {}; Jdl: {}; Jtl: {}'.format(len(Jsl), len(Jdl), len(Jtl)))

def main():
    # command-line argument parser
    parser = argparse.ArgumentParser(description='Identification of synthetic lethal reactions of a genome-wide metabolic network model')
    parser.add_argument('model',nargs='?',type=str,help='model .xml file')
    parser.add_argument('--cutoff',nargs='?',default=0.01,type=int,help='cut-off value')
    parser.add_argument('--order',nargs='?',default=2,type=int,help='order')
    parser.add_argument('--eliList',nargs='?',type=str,help='.xml file of elimination reactions')
    parser.add_argument('--atpm',nargs='?',default='ATPM',type=str,help='ID of ATP maintenance reaction')
    parser.add_argument('--solver',nargs='?',default='glpk',type=str,help='solver used')

    args = parser.parse_args()

    # Logger initiation
    logging.getLogger().setLevel(logging.INFO)

    # FastSL
    fastSL(model=args.model,cutoff=args.cutoff,order=args.order,eliList=args.eliList,atpm=args.atpm,solver=args.solver)

if __name__ == '__main__':
    main()
