# -*- coding: utf-8 -*-

import csv
from pathlib import Path
import logging

from lxml import etree


logger = logging.getLogger(__name__)


def generate_dir(model_name, of_type):
    '''Generate the required directory for output'''
    dir_path = Path('Results/{}/{}'.format(model_name, of_type))
    dir_path.mkdir(parents=True, exist_ok=True)

    logger.info('%s generated.', str(dir_path.resolve()))
    return dir_path


def write_file(in_dir, model_name, of_type, lethality_stage, data):
    file_path = in_dir / '{}_{}_lethal_{}.csv'.format(model_name,
                                                      lethality_stage,
                                                      of_type)

    file_path.touch()
    logger.info('%s generated.', str(file_path.resolve()))

    with file_path.open('w') as file:
        if lethality_stage == 'single':
            csv.writer(file).writerow(data)
        elif lethality_stage == 'double' or lethality_stage == 'triple':
            csv.writer(file).writerows(data)


def handle_elilist(model, atpm, elilist):
    if not elilist:
        # error handling for adding ATP maintenance reaction
        # in the elimination list if no elilist provided
        try:
            elilist_data = model.reactions.get_by_id(atpm).id
            logger.info('%s found in model and added to elimination list.',
                        atpm)
        except KeyError:
            logger.info('%s not found in the model.', atpm)

    elif elilist:
        # elimination list parsing
        elilist_tree = etree.parse(elilist)
        elilist_data = [data.text for data in
                        elilist_tree.iter(tag='reaction-id')]

        return elilist_data
