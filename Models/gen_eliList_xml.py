# -*- coding: utf-8 -*-

import argparse
import logging
import os

from cobra.io import read_sbml_model
from lxml import etree


def generate_elilist_xml_from_model(model):
    '''
    Generates .xml file for the elimination list used for analysis of
    synthetic lethals for a given models
    '''

    # Exchange reactions -> Indices -> Reaction IDs
    model = read_sbml_model(model)
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

    logging.info('{}_elimination_list.xml generated in {}'
                 .format(model, os.getcwd()))


def main():
    '''
    Takes command line arguments and generates the .xml file
    '''
    # command-line argument parser
    parser = argparse.ArgumentParser(description='Generation of reaction\
                                                  elimination list .xml file\
                                                  for a genome-wide\
                                                  metabolic network model')
    parser.add_argument('model',
                        nargs='?',
                        type=str,
                        help='model .xml file')
    args = parser.parse_args()

    # Logger initiation
    logging.getLogger().setLevel(logging.INFO)

    generate_elilist_xml_from_model(args.model)


if __name__ == '__main__':
    main()
