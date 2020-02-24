#!/usr/bin/env python3
import re
from os import path, mkdir
import logging
from re import search
from subprocess import check_output

from yaml import safe_load

from src.load_barcode import load_barcode

import src.preprocessing as s0
import src.fastq2fasta as s1
import src.remove_primer as s2
import src.remove_vector as s3
import src.velvet_assembly as s4
import src.phrap_assembly as s5


def init_logger(work_dir):
    if not path.exists(work_dir):
        mkdir(work_dir)
    logging.basicConfig(
        filename=path.join(work_dir, 'parfums.log'),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


if __name__ == '__main__':
    with open('examples/configs.yml') as yaml_file:
        config = safe_load(yaml_file)
    work_dir = config['work_dir']
    primer_file = config['primer_file']
    vector_file = config['vector_file']
    fw_file = config['FW_file']
    rc_file = config['RC_file']
    bc_file = config['BC_file']

    init_logger(work_dir)
    barcodes, length = load_barcode(bc_file)
    idents = list(barcodes.values())

    # s0.main(fw_file, rc_file, work_dir, barcodes, length, mismatch=False)
    # s1.main(work_dir, idents)
    # s2.clean_seq(work_dir, primer_file, idents, 200000)
    # s3.main(work_dir, vector_file, idents, 200000)
    # s4.main(work_dir, idents)
    # s5.main(work_dir, idents)
