#!/usr/bin/env python3
import logging
from os import path


def merge_file(work_dir, split_files, suffix):
    logging.info('Merging clean fasta output files')
    result = dict()
    for (ident, files) in split_files.items():
        out_file = path.join(work_dir, ident, '{}.{}'.format(ident, suffix))
        result[ident] = out_file
        with open(out_file, 'w') as out:
            for file in files:
                with open(file) as split_cm:
                    for line in split_cm:
                        out.write(line)
        logging.info('{} is formed'.format(out_file))
    return result
