#!/usr/bin/env python3
import logging
from sys import path


def merge_file(work_dir, split_files, suffix):
    logging.info('Merging clean fasta output files')
    merged_file = dict()
    for (ident, files) in split_files.items():
        out_file = path.join(work_dir, ident, '{}.{}'.format(ident, suffix))
        merged_file[ident] = out_file
        with open(out_file, 'w') as out:
            for file in files:
                with open(file) as split_file:
                    for line in split_file:
                        out.write(line)
        logging.info('{} is merged'.format(out_file))
    return merged_file
