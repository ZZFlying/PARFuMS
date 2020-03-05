#!/usr/bin/env python3
import sys
import logging
import traceback
from os import path, mkdir

from src.load_barcode import load_barcode

from src.preprocessing import main as preprocessing
from src.fastq2fasta import main as fastq2fasta
from src.remove_primer import main as remove_primer
from src.remove_vector import main as remove_vector
from src.velvet_assembly import main as velvet_assembly
from src.phrap_assembly import main as phrap_assembly


def init_logger(work_dir):
    logging.basicConfig(
        filename=path.join(work_dir, 'parfums.log'),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s %(module)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_pre_file(work_dir, idents):
    fastq_dir = path.join(work_dir, 'PreprocessedFiles')
    for ident in idents:
        fw_file = path.join(fastq_dir, '{}_FW.fq'.format(ident))
        rc_file = path.join(fastq_dir, '{}_RC.fq'.format(ident))
        if path.exists(fw_file) and path.exists(rc_file):
            continue
        idents.remove(ident)
    return idents


def check_output_files(workdir, idents, step, suffix):
    not_pass_list = list()
    for ident in idents:
        out_file = path.join(workdir, ident, '{}.{}'.format(ident, suffix))
        if path.exists(out_file) and path.getsize(out_file):
            continue
        logging.info('{} is empty/not found. Halting its processing here'.format(out_file))
        not_pass_list.append(ident)
        idents.remove(ident)
    if not_pass_list:
        not_pass_seq = ''.join([i for i in not_pass_list])
        logging.info('IDs not included at Step{}: {}'.format(step, not_pass_seq))
    return idents

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


def read_config(config_file):
    with open(config_file) as config_in:
        contig = dict()
        for line in config_in:
            if line.startswith('#'):
                continue
            line = line.strip()
            array = line.split(':')
            key = array[0].strip()
            value = array[1].strip()
            if value.isdigit():
                value = int(value)
            elif value.upper() == 'TRUE':
                value = True
            elif value.upper() == 'FALSE':
                value = False
            contig[key] = value
        return contig


def parse_step(step):
    step = step.strip()
    step_list = list()
    if ':' in step:
        step = step.split(':')
        start = int(step[0])
        stop = int(step[1])
        step_list = range(start, stop + 1)
    elif int(step) == 0:
        step_list = range(1, 6)
    else:
        step_list = [int(step)]
    return step_list

if __name__ == '__main__':
    # config_file = sys.argv[1]
    # step = sys.argv[2]
    # config = read_config(config_file)
    # step_list = parse_step(step)

    config = read_config('examples/configs.yml')
    step_list = parse_step('5')

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
    try:
        for step in step_list:
            if step == 1:
                preprocessing(fw_file, rc_file, work_dir, barcodes, length, mismatch)
                idents = check_pre_file(work_dir, idents)
                fastq2fasta(work_dir, idents, maxsize)
                logging.info('STEP 1: PREPROCESSING COMPLETED')
            elif step == 2:
                idents = check_output_files(work_dir, idents, step, suffix='fasta')
                remove_primer(work_dir, primer_file, idents)
                logging.info('STEP 2: PRIMER SEQUENCES REMOVED')
            elif step == 3:
                idents = check_output_files(work_dir, idents, step, suffix='clean.fasta')
                remove_vector(work_dir, vector_file, idents)
                logging.info('STEP 3: VECTOR SEQUENCES REMOVED')
            elif step == 4:
                idents = check_output_files(work_dir, idents, step, suffix='noVector.fasta')
                velvet_assembly(work_dir, idents)
                logging.info('STEP 4: VELVET ASSEMBLY COMPLETED')
            elif step == 5:
                idents = check_output_files(work_dir, idents, step, suffix='ForPhrap1.fasta')
                phrap_assembly(work_dir, idents)
                logging.info('STEP 5: PHRAP ASSEMBLY COMPLETED')
            else:
                print('Step Count doesn\'t make sense in current program')
                logging.info('Step Count doesn\'t make sense in current program')
                exit(1)
    except KeyboardInterrupt:
        logging.info('INTERRUPT THE PARFUMS PIPELINE')
    except Exception as e:
        traceback.print_exc()
        logging.info('EXIT WITH ERROR')
        logging.info(traceback.format_exc())
    logging.info('END OF PARFUMS PIPELINE')

