#!/usr/bin/env python3
import logging
from os import path, mkdir
from yaml import safe_load

from src.load_barcode import load_barcode

from src.preprocessing import main as preprocessing
from src.fastq2fasta import main as fastq2fasta
from src.remove_primer import main as remove_primer
from src.velvet_assembly import main as velvet_assembly
import src.phrap_assembly


def init_logger(work_dir):
    logging.basicConfig(
        filename=path.join(work_dir, 'parfums.log'),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s %(module) - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_output_files(workdir, idents, suffix):
    uncheck = list()
    for ident in idents:
        out_file = path.join(workdir, ident, ident + suffix)
        if not path.exists(out_file) or not path.getsize(out_file):
            logging.info('{} is empty/not found. Halting its processing here'.format(out_file))
            uncheck.append(ident)
    return uncheck

def make_output_dir(work_dir):
    preprocessed_files = path.join(work_dir, 'PreprocessedFiles')
    temp = path.join(work_dir, 'temp')
    sub_script = path.join(temp, 'qsub_scripts')
    if not path.exists(work_dir):
        mkdir(work_dir)
        mkdir(preprocessed_files)
        mkdir(temp)
        mkdir(sub_script)
    elif not path.exists(preprocessed_files):
        mkdir(preprocessed_files)
    elif not path.exists(temp):
        mkdir(temp)
    elif not path.exists(sub_script):
        mkdir(sub_script)

if __name__ == '__main__':
    with open('examples/configs.yml') as yaml_file:
        config = safe_load(yaml_file)
    work_dir = config['work_dir']
    primer_file = config['primer_file']
    vector_file = config['vector_file']
    fw_file = config['FW_file']
    rc_file = config['RC_file']
    bc_file = config['BC_file']
    mismatch = config['mismatch']
    maxsize = config['maxsize']
    make_output_dir(work_dir)
    init_logger(work_dir)
    barcodes, idents, length = load_barcode(bc_file)

    # preprocessing(fw_file, rc_file, work_dir, barcodes, length, mismatch)
    # fastq2fasta(work_dir, idents, maxsize)
    # remove_primer(work_dir, primer_file, vector_file, idents)
    velvet_assembly(work_dir, idents)
