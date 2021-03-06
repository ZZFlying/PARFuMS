#!/usr/bin/env python3

import logging
from os import path
from shutil import rmtree
from tempfile import mkdtemp

from util.merge_split import merge_file
from util.compress import compress_multi
from util.singleton_config import Config
from util.scriptReator import cross_match_script

from sub.parfums_subs import submit_array, read_fasta


def main(work_dir, primer_file, idents, is_gzip, maxsize=200000):
    logging.info('STEP 2: REMOVING PRIMER SEQUENCES')
    tempdir = mkdtemp(dir=path.join(work_dir, 'temp'), prefix='primer_')
    script = 'CleanPrimer.py'
    # 将大文件分割为maxsize条序列的小文件，多进程同时进行任务
    suffix = 'fasta.gz' if is_gzip else 'fasta'
    split_file = read_fasta(work_dir, tempdir, idents, maxsize, suffix)
    script, result = cross_match_script(split_file, primer_file, tempdir, script, suffix='clean')
    submit_array(script, 'primer_run', tempdir)
    # 合并清除序列后的分割文件
    files_list = merge_file(work_dir, result, suffix='clean.fasta')
    compress_multi(files_list.values())
    if Config()['auto_del']:
        rmtree(tempdir)
