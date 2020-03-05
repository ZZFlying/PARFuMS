#!/usr/bin/env python3

import sys
import logging
from os import path
from collections import defaultdict

from sub.parfums_subs import submit_array, create_tempdir, read_fasta


def make_script(split_file, map_file, tempdir, script, suffix):
    cm_files = defaultdict(list)
    cm_script = 'CM_run_script.sh'
    cm_script = path.join(tempdir, cm_script)
    sub_folder = path.join(sys.path[0], 'sub')
    script = path.join(sub_folder, script)

    with open(cm_script, 'w') as output:
        for (ident, files) in split_file.items():
            for file in files:
                cmd = 'cross_match {} {}'.format(file, map_file)
                cmd += ' -gap1_only -minmatch 6 -minscore 10 -gap_init -3'
                # cmd += ' | perl {} {}'.format(script, file)
                cmd += ' | python3 {} {}'.format(script, file)
                cmd += '\n'
                output.write(cmd)
                cm_file = '{}.{}'.format(file, suffix)
                cm_files[ident].append(cm_file)
    return cm_script, cm_files


def merge_cm_file(work_dir, cm_files, suffix):
    logging.info('Merging clean fasta output files')
    merge_file = dict()
    for (ident, files) in cm_files.items():
        out_file = path.join(work_dir, ident, '{}.{}'.format(ident, suffix))
        merge_file[ident] = out_file
        with open(out_file, 'w') as out:
            for file in files:
                with open(file) as split_cm:
                    for line in split_cm:
                        out.write(line)
        logging.info('{} is formed'.format(out_file))
    return merge_file


def clean_primer_seq(work_dir, primer_file, idents, maxsize):
    logging.info('STEP 2: REMOVING PRIMER SEQUENCES')
    directory = create_tempdir(work_dir, prefix='primer_CMRun_')
    tempdir = directory.name
    # script = 'CleanAdapter-Illumina_PE_mod.pl'
    script = 'CleanAdapter-Illumina_PE_mod.py'
    split_file = read_fasta(work_dir, tempdir, idents, maxsize, suffix='fasta')
    cm_script, cm_files = make_script(split_file, primer_file, tempdir, script, suffix='clean')
    submit_array(cm_script, 'primer_CM_run', tempdir)
    merge_cm_file(work_dir, cm_files, suffix='clean.fasta')


def main(work_dir, primer_file, idents, maxsize=200000):
    clean_primer_seq(work_dir, primer_file, idents, maxsize)
