# -*- coding: utf-8 -*-

import logging
import sys

import cobra

from lxml import etree

from fast_sl.util import generate_dir, write_file, cli_parser


def fast_sl(model_path, cutoff, order, atpm, solver, parallel, elilist=None):
    '''Identify synthetic lethal reactions.'''
    # model parsing
    model = cobra.io.read_sbml_model(model_path)

    # directory creation for output
    results_dir_path = generate_dir(model.id, "reactions")

    # elilist data processing
    if not elilist:
        # error handling for adding ATP maintenance reaction
        # in the elimination list if no elilist provided
        try:
            elilist_data = [model.reactions.get_by_id(atpm).id]
            logging.info('%s found in model and added to elimination list.',
                         atpm)
        except KeyError:
            logging.critical(('%s not found in model. Enter a valid ATP '
                              'maintenance reaction ID or provide an '
                              'elimination list.'), atpm)
            sys.exit()

    elif elilist:
        # elimination list parsing
        elilist_tree = etree.parse(elilist)
        elilist_data = [data.text for data in
                        elilist_tree.iter(tag='reaction-id')]

    # core handling
    if parallel == 0:
        from fast_sl.single_core import single_sl, double_sl
        # Order 1
        if order == 1:
            jsl = single_sl(model,
                            cutoff,
                            elilist_data,
                            solver)
            write_file(results_dir_path,
                       model.id,
                       "reactions",
                       "single",
                       jsl)
            logging.info('========== %s ==========', model.id)
            logging.info('Jsl:%s', len(jsl))
        # Order 2
        elif order == 2:
            jsl, jdl = double_sl(model,
                                 cutoff,
                                 elilist_data,
                                 solver)
            write_file(results_dir_path,
                       model.id,
                       "reactions",
                       "single",
                       jsl)
            write_file(results_dir_path,
                       model.id,
                       "reactions",
                       "double",
                       jdl)
            logging.info('========== %s ==========', model.id)
            logging.info('Jsl:%s\tJdl:%s', len(jsl), len(jdl))
        # Order 3
        # elif order == 3:
        #     jsl, jdl, jtl = triple_sl(model,
        #                               cutoff,
        #                               elilist_data,
        #                               solver)
        #     write_file(results_dir_path,
        #                model.id,
        #                "reactions",
        #                "single",
        #                jsl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "reactions",
        #                "double",
        #                jdl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "reactions",
        #                "triple",
        #                jtl)
        #     logging.info('Jsl:%s\tJdl:%s\tJtl:%s',
        #                  len(jsl),
        #                  len(jdl),
        #                  len(jtl))
    elif parallel == 1:
        from fast_sl.multi_core import (
                                        parallel_single_sl,
                                        parallel_double_sl)
        # Parallel order 1
        if order == 1:
            jsl = parallel_single_sl(model,
                                     cutoff,
                                     elilist_data,
                                     solver)
            write_file(results_dir_path,
                       model.id,
                       "reactions",
                       "single",
                       jsl)
            logging.info('========== %s ==========', model.id)
            logging.info('Jsl:%s', len(jsl))
        # Parallel order 2
        elif order == 2:
            jsl, jdl = parallel_double_sl(model,
                                          cutoff,
                                          elilist_data,
                                          solver)
            write_file(results_dir_path,
                       model.id,
                       "reactions",
                       "single",
                       jsl)
            write_file(results_dir_path,
                       model.id,
                       "reactions",
                       "double",
                       jdl)
            logging.info('========== %s ==========', model.id)
            logging.info('Jsl:%s\tJdl:%s',
                         len(jsl),
                         len(jdl))
        # Parallel order 3
        # elif order == 3:
        #     jsl, jdl, jtl = parallel_triple_sl(model,
        #                                        cutoff,
        #                                        elilist_data,
        #                                        solver)
        #     write_file(results_dir_path,
        #                model.id,
        #                "reactions",
        #                "single",
        #                jsl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "reactions",
        #                "double",
        #                jdl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "reactions",
        #                "triple",
        #                jtl)
        #     logging.info('jsl:%s\tjdl:%s\tjtl:%s',
        #                 len(jsl),
        #                 len(jdl),
        #                 len(jtl))


def fast_sl_genes(model_path, cutoff, order, solver, parallel):
    '''Identify synthetic lethal genes.'''
    # model parsing
    model = cobra.io.read_sbml_model(model_path)

    # directory creation for output
    results_dir_path = generate_dir(model.id, "genes")

    # core handling
    if parallel == 0:
        from fast_sl.single_core import (
                                         single_sl_genes,
                                         double_sl_genes)
        # Order 1
        if order == 1:
            jsl = single_sl_genes(model,
                                  cutoff,
                                  solver)
            write_file(results_dir_path,
                       model.id,
                       "genes",
                       "single",
                       jsl)
            logging.info('========== %s ==========', model.id)
            logging.info('jsl genes:%s', len(jsl))
        # Order 2
        elif order == 2:
            jsl, jdl = double_sl_genes(model,
                                       cutoff,
                                       solver)
            write_file(results_dir_path,
                       model.id,
                       "genes",
                       "single",
                       jsl)
            write_file(results_dir_path,
                       model.id,
                       "genes",
                       "double",
                       jdl)
            logging.info('========== %s ==========', model.id)
            logging.info('Jsl genes:%s\tJdl genes:%s', len(jsl), len(jdl))
        # Order 3
        # elif order == 3:
        #     jsl, jdl, jtl = triple_sl_genes(model,
        #                                     cutoff,
        #                                     elilist_data,
        #                                     solver)
        #     write_file(results_dir_path,
        #                model.id,
        #                "genes",
        #                "single",
        #                jsl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "genes",
        #                "double",
        #                jdl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "genes",
        #                "triple",
        #                jtl)
        #     logging.info('jsl genes:%s\tjdl genes:%s\tjtl genes:%s',
        #                 len(jsl),
        #                 len(jdl),
        #                 len(jtl))
    elif parallel == 1:
        from fast_sl.multi_core import (
                                        parallel_single_sl_genes,
                                        parallel_double_sl_genes)
        # Parallel order 1
        if order == 1:
            jsl = parallel_single_sl_genes(model,
                                           cutoff,
                                           solver)
            write_file(results_dir_path,
                       model.id,
                       "genes",
                       "single",
                       jsl)
            logging.info('========== %s ==========', model.id)
            logging.info('Jsl genes:%s', len(jsl))
        # Parallel order 2
        elif order == 2:
            jsl, jdl = parallel_double_sl_genes(model,
                                                cutoff,
                                                solver)
            write_file(results_dir_path,
                       model.id,
                       "genes",
                       "single",
                       jsl)
            write_file(results_dir_path,
                       model.id,
                       "genes",
                       "double",
                       jdl)
            logging.info('========== %s ==========', model.id)
            logging.info('Jsl genes:%s\tJdl genes:%s',
                         len(jsl),
                         len(jdl))
        # Parallel order 3
        # elif order == 3:
        #     jsl, jdl, jtl = parallel_triple_sl_genes(model,
        #                                              cutoff,
        #                                              solver)
        #     write_file(results_dir_path,
        #                model.id,
        #                "genes",
        #                "single",
        #                jsl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "genes",
        #                "double",
        #                jdl)
        #     write_file(results_dir_path,
        #                model.id,
        #                "genes",
        #                "triple",
        #                jtl)
        #     logging.info('jsl genes:%s\tjdl genes:%s\tjtl genes:%s',
        #                 len(jsl),
        #                 len(jdl),
        #                 len(jtl))


def main():
    '''Take parameters and process.'''
    # command-line arguments
    args = cli_parser()

    # logging
    logging.basicConfig(filename='fast-sl.log',
                        format='%(asctime)s : %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)

    # elilist generation handling
    if args.gen_elilist is True:
        from fast_sl.util import generate_elilist_xml_from_model
        generate_elilist_xml_from_model(args.model)
    else:
        if args.genes == 0:
            # Fast-sl
            fast_sl(args.model,
                    args.cutoff,
                    args.order,
                    args.atpm,
                    args.solver,
                    args.parallel,
                    args.elilist)
        else:
            # Fast-sl genes
            fast_sl_genes(args.model,
                          args.cutoff,
                          args.order,
                          args.solver,
                          args.parallel)


if __name__ == '__main__':
    main()
