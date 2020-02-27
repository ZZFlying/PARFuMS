#!/usr/bin/env python3

import sys, logging

from os import path, mkdir
from shutil import copy
from collections import defaultdict
from subprocess import check_output

from sub.parfums_subs import submit_array, split_fasta_files, create_tempdir


def read_fasta(work_dir, tempdir, idents, maxsize, suffix):
    split_file = dict()
    process_file = dict()
    split_dir = path.join(tempdir, 'splitFastaFiles')
    if not path.exists(split_dir):
        mkdir(split_dir)
    for ident in idents:
        fasta_filename = ident + suffix
        fasta_file = path.join(work_dir, ident, fasta_filename)
        if path.exists(fasta_file):
            cmd = "grep '^>' {} | wc -l".format(fasta_file)
            count = check_output(cmd, shell=True).decode('utf-8')
            count = int(count)
            if count > maxsize:
                process_file[ident] = fasta_file
            else:
                copy(fasta_file, split_dir)
                split_file[ident] = path.join(split_dir, fasta_filename)
    process_file = split_fasta_files(process_file, split_dir, maxsize)
    split_file.update(process_file)
    return split_file


def make_crossmatch_script(split_file, vector_file, tempdir, script, suffix):
    cm_files = defaultdict(list)
    cm_script = 'CM_run_script.sh'
    cm_script = path.join(tempdir, cm_script)
    sub_folder = path.join(sys.path[0], 'sub')
    script = path.join(sub_folder, script)
    logging.info('Making jobscript file')
    with open(cm_script, 'w') as output:
        for (ident, files) in split_file.items():
            for file in files:
                cmd = 'cross_match {} {}'.format(file, vector_file) \
                      + ' -gap1_only -minmatch 6 -minscore 10 -gap_init -3' \
                      + ' | ' \
                      + 'perl {} {}\n'.format(script, file)
                output.write(cmd)
                cm_file = file + suffix
                cm_files[ident].append(cm_file)
    return cm_script, cm_files


def merge_cm_file(work_dir, cm_files, suffix):
    logging.info('Merging clean fasta output files')
    outfile = ''
    for (ident, files) in cm_files.items():
        outfile = path.join(work_dir, ident, ident + suffix)
        with open(outfile, 'w') as out:
            for file in files:
                with open(file) as clean:
                    for line in clean:
                        out.write(line)
        logging.info('{} is formed'.format(outfile))
    return outfile


def main(work_dir, vector_file, idents, maxsize):
    logging.info('STEP 3: REMOVING VECTOR SEQUENCES')
    directory = create_tempdir(work_dir, prefix='vector_CMRun_')
    tempdir = directory.name
    split_file = read_fasta(work_dir, tempdir, idents, maxsize, suffix='.clean.fasta')
    cm_script, cm_files = make_crossmatch_script(split_file, vector_file, tempdir,
                                                 script='RemoveVector_Gautam.pl',
                                                 suffix='.noVector')
    submit_array(cm_script, 'vector_CM_run', tempdir)
    merge_cm_file(work_dir, cm_files, suffix='.noVector.fasta')