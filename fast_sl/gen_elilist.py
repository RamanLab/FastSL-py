# -*- coding: utf-8 -*-

import logging
import os

from cobra.io import read_sbml_model
from lxml import etree


logger = logging.getLogger(__name__)


def generate_elilist_xml_from_model(model):
    ''' Generates .xml file for the elimination list used for analysis of
    synthetic lethals for a given models.'''
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
    tree.write('Models/{}_elimination_list.xml'.format(model),
               encoding='utf-8',
               xml_declaration=True,
               pretty_print=True)

    logger.info('%s_elimination_list.xml generated in Models',
                model)
