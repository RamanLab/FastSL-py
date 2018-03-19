# -*- coding: utf-8 -*-

import argparse


def cli_parser():
    '''Take command line arguments.'''
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
    parser.add_argument('--parallel',
                        nargs='?',
                        default=0,
                        type=int,
                        help='0: single-core, 1: multi-core')
    parser.add_argument('--gen_elilist',
                        nargs='?',
                        default=False,
                        type=bool,
                        help='model .xml file')
    args = parser.parse_args()

    return args
