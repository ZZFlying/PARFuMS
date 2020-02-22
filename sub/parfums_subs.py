#!/usr/bin/env python3

import sys, logging

from re import search
from os import path, remove, mkdir
from glob import glob
from math import ceil
from shutil import copy
from tempfile import TemporaryDirectory
from time import sleep
from subprocess import check_output, CalledProcessError


def create_tempdir(work_dir, prefix):
    tempdir = path.join(work_dir, 'temp')
    if not path.exists(tempdir):
        mkdir(tempdir)
    directory = TemporaryDirectory(prefix=prefix, dir=tempdir)
    return directory


def read_fasta(work_dir, tempdir, idents, maxsize, suffix):
    split_file = dict()
    process_file = dict()
    split_dir = path.join(tempdir, 'splitFastaFiles')
    if not path.exists(split_dir):
        mkdir(split_dir)
    logging.info('Value of MaxSize inside make vector run: {}'.format(maxsize))
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
                split_file[ident] = [path.join(split_dir, fasta_filename)]
    process_file = split_fasta_files(process_file, split_dir, maxsize)
    split_file.update(process_file)
    return split_file


def split_fasta_files(process_file, out_dir, maxsize):
    results = dict()
    for (ident, file) in process_file.items():
        split_list = list()
        with open(file) as f:
            count = 0
            file_part = 1
            split_name = '{}-part{}.fasta'
            split_file = path.join(out_dir, split_name.format(ident, file_part))
            output = open(split_file, 'w')
            line = f.readline()
            while line:
                if line.startswith('>'):
                    if count < maxsize:
                        output.writelines([line, f.readline()])
                        count += 1
                    else:
                        output.close()
                        split_list.append(split_file)

                        count = 0
                        file_part += 1
                        split_file = path.join(out_dir, split_name.format(ident, file_part))
                        output = open(split_file, 'w')
                        output.writelines([line, f.readline()])
                        count += 1
                line = f.readline()
            output.close()
            split_list.append(split_file)
        results[ident] = split_list
    return results


def submit_array(job_script, job_name, work_dir):
    with open(job_script) as file:
        job_count = len(file.readlines())
    step_size = ceil(job_count / 1000)
    sub_folder = path.join(sys.path[0], 'sub')
    slurm_script = path.join(sub_folder, 'submit_job_array_slurm.sh')
    slurm_cmd = 'sbatch --job-name={} --chdir={} --mem=2000'.format(job_name, work_dir)
    slurm_cmd = str.join(' ', [slurm_script, job_script, str(step_size), '|', slurm_cmd])
    logging.debug('cmd =>' + slurm_cmd)

    job_return = check_output(slurm_cmd, shell=True).decode('utf-8')
    logging.debug('job_return =>' + job_return.strip('\n'))
    job_id = search('Submitted batch job (\d+)', job_return).group(1)

    completed_file = '{}_{}_check.txt'.format(job_name, job_id)
    completed_file = path.join(work_dir, completed_file)
    complete_script = path.join(sub_folder, 'complete.sh')
    if path.exists(completed_file):
        remove(completed_file)

    slurm_cmd = "sbatch --dependency=singleton " \
                + "--chdir={} ".format(work_dir) \
                + "--job-name={} {} {} | ".format(job_name, complete_script, completed_file) \
                + "awk '{print $NF}'"
    logging.debug('cmd =>' + slurm_cmd)

    job_return = check_output(slurm_cmd, shell=True).decode('utf-8')
    logging.debug('job_return =>' + job_return.strip('\n'))

    wait_job_finished(job_id, work_dir, completed_file, job_count)


def wait_job_finished(job_id, work_dir, completed_file, job_count):
    logging.info('Waiting for job: {} to finish!'.format(job_id))
    suc_count = 0
    err_count = 0
    while True:
        if path.exists(completed_file):
            logging.info('job {} finished...Will now check for errors'.format(job_id))
            slurm_out = path.join(work_dir, 'slurm-{}_*'.format(job_id))
            files = glob(slurm_out)
            for file in files:
                suc_cmd = "grep -h -c 'Well done!! Job finished' {}".format(file)
                try:
                    check_output(suc_cmd, shell=True)
                    suc_count += 1
                except CalledProcessError:
                    err_count += 1
            logging.info('success_cnt => {}'.format(suc_count))
            logging.info('error_cnt => {}'.format(err_count))
            break
        sleep(5)
    if suc_count == job_count and err_count == 0:
        logging.info('job {} completed successfully'.format(job_id))
    else:
        logging.info('Error occurred while running the job {}'.format(job_id))
