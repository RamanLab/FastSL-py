#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import argparse, logging, os
from cobra.io import read_sbml_model
from lxml import etree

def generate_elilist_xml_from_model(model):
    '''
    Generates .xml file for the elimination list used for analysis of synthetic lethals for a given model
    '''
    model = read_sbml_model(model)
    exchange_reactions_objects = model.exchanges # cobra.core.Reaction objects
    exchange_reactions_index = [model.reactions.index(reaction_id) for reaction_id in exchange_reactions_objects]
    exchange_reactions_list = [model.reactions[reaction_idx].id for reaction_idx in exchange_reactions_index]

    # xml generation
    root = etree.Element('elimination-list')
    for reaction_id in exchange_reactions_list:
        child = etree.SubElement(root, 'reaction-id')
        child.text = reaction_id
    tree = etree.ElementTree(root)
    tree.write('{}_elimination_list.xml'.format(model),encoding='utf-8',xml_declaration=True,pretty_print=True)

    logging.info('{}_elimination_list.xml generated in {}'.format(model, os.getcwd()))

def main():
    # command-line argument parser
    parser = argparse.ArgumentParser(description='Generation of elimination list .xml file for a genome-wide metabolic network model')
    parser.add_argument('model',nargs='?',type=str,help='model .xml file')
    args = parser.parse_args()

    # Logger initiation
    logging.getLogger().setLevel(logging.INFO)

    generate_elilist_xml_from_model(model=args.model)

if __name__ == '__main__':
    main()