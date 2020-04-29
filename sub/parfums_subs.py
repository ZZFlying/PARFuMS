#!/usr/bin/env python3

import sys
import logging
from multiprocessing import Pool, Manager

from re import search
from os import path, remove, mkdir
from glob import glob
from math import ceil
from shutil import copy
from tempfile import TemporaryDirectory
from time import sleep

import gzip

from util.singleton_config import Config
from subprocess import check_output, CalledProcessError, Popen


def create_tempdir(work_dir, prefix):
    tempdir = path.join(work_dir, 'temp')
    if not path.exists(tempdir):
        mkdir(tempdir)
    directory = TemporaryDirectory(prefix=prefix, dir=tempdir)
    logging.info('Created temp folder: {}'.format(tempdir))
    return directory


def read_fasta(work_dir, tempdir, idents, maxsize, suffix):
    split_file = dict()
    process_file = dict()
    split_dir = path.join(tempdir, 'splitFastaFiles')
    if not path.exists(split_dir):
        mkdir(split_dir)
    logging.info('Value of MaxSize inside make vector run: {}'.format(maxsize))
    for ident in idents:
        fasta_filename = '{}.{}'.format(ident, suffix)
        fasta_file = path.join(work_dir, ident, fasta_filename)
        if path.exists(fasta_file):
            cmd_template = 'zcat {} | grep "^>" | wc -l' if suffix.split('.')[-1] == 'gz' else "grep '^>' {} | wc -l"
            cmd = cmd_template.format(fasta_file)
            count = int(check_output(cmd, shell=True))
            if count > maxsize:
                process_file[ident] = fasta_file
            else:
                copy(fasta_file, split_dir)
                split_file[ident] = [path.join(split_dir, fasta_filename)]
    logging.debug('process list:' + process_file.__str__())
    process_file = split_fasta_files(process_file, split_dir, maxsize)
    split_file.update(process_file)
    return split_file


def split_fasta_files(process_file, out_dir, maxsize):
    results = dict()
    for (ident, file) in process_file.items():
        split_list = list()
        _open = gzip.open if file.split('.')[-1] == 'gz' else open
        with _open(file, 'rt') as f:
            count = 0
            file_part = 1
            split_name = '{}-part{}.fasta'
            split_file = path.join(out_dir, split_name.format(ident, file_part))
            output = open(split_file, 'wt')
            line = f.readline()
            while line:
                if line.startswith('>'):
                    if count < maxsize:
                        output.write(line)
                        count += 1
                    else:
                        output.close()
                        split_list.append(split_file)

                        count = 0
                        file_part += 1
                        split_file = path.join(out_dir, split_name.format(ident, file_part))
                        output = open(split_file, 'w')
                        output.write(line)
                        count += 1
                else:
                    output.write(line)
                line = f.readline()
            output.close()
            split_list.append(split_file)
        results[ident] = split_list
    return results


def executor_foo(task_queue, done_queue, Id, work_dir):
    try:
        out_file = 'mission-{}.out'.format(Id)
        out_file = path.join(work_dir, out_file)
        out_file = open(out_file, 'w')
        mission = task_queue.get()
        Popen(mission, shell=True, stdout=out_file, stderr=out_file).wait()
        done_queue.put(True)
    except CalledProcessError as e:
        logging.ERROR(e.cmd)
        done_queue.put(False)


def submit_array(job_script, job_name, work_dir):
    with open(job_script) as file:
        job_list = list()
        for line in file:
            line = line.strip()
            job_list.append(line)
    task_queue = Manager().Queue()
    done_queue = Manager().Queue()
    threads = Config()['thread']
    logging.info('Working with {} threads'.format(threads))
    with Pool(processes=threads) as executor:
        logging.info('Waiting for job: {} to finish!'.format(job_name))
        # for line in job_list:
        for i in range(len(job_list)):
            task_queue.put(job_list[i])
            executor.apply_async(executor_foo, (task_queue, done_queue, i, work_dir))
        finish_jobs = 0
        executor.close()
        executor.join()
        logging.info('job {} finished...Will now check for errors'.format(job_name))
        for i in range(len(job_list)):
            result = done_queue.get()
            if result:
                finish_jobs += 1
        logging.info('success_cnt => {}'.format(finish_jobs))
        logging.info('error_cnt => {}'.format(len(job_list) - finish_jobs))


def submit_array_old(job_script, job_name, work_dir):
    with open(job_script) as file:
        job_count = len(file.readlines())
    step_size = ceil(job_count / 1000)
    sub_folder = path.join(sys.path[0], 'sub')
    slurm_script = path.join(sub_folder, 'submit_job_array_slurm.sh')
    slurm_cmd = 'sbatch -J {} -D {} --mem=2000'.format(job_name, work_dir)
    slurm_cmd = str.join(' ', [slurm_script, job_script, str(step_size), '|', slurm_cmd])
    logging.debug('cmd => ' + slurm_cmd)

    job_return = check_output(slurm_cmd, shell=True).decode('utf-8')
    logging.debug('job_return => ' + job_return.strip())
    job_id = search('Submitted batch job (\d+)', job_return).group(1)

    completed_file = '{}_{}_check.txt'.format(job_name, job_id)
    completed_file = path.join(work_dir, completed_file)
    complete_script = path.join(sub_folder, 'complete.sh')
    if path.exists(completed_file):
        remove(completed_file)

    slurm_cmd = "sbatch --dependency=singleton " \
                + "-D {} ".format(work_dir) \
                + "-J {} {} {} | ".format(job_name, complete_script, completed_file) \
                + "awk '{print $NF}'"
    logging.debug('cmd => ' + slurm_cmd)

    job_return = check_output(slurm_cmd, shell=True).decode('utf-8')
    logging.debug('job_return => ' + job_return.strip())

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
