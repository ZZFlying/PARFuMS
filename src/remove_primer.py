#!/usr/bin/env python3

import sys
import logging
from os import path
from collections import defaultdict
from shutil import rmtree
from tempfile import mkdtemp

from util.compress import compress_multi
from util.merge_split import merge_file
from sub.parfums_subs import submit_array, read_fasta
from util.singleton_config import Config


def make_script(split_file, map_file, tempdir, script, suffix):
    cm_files = defaultdict(list)
    cm_script = 'CM_run_script.sh'
    cm_script = path.join(tempdir, cm_script)
    sub_folder = path.join(sys.path[0], 'sub')
    script = path.join(sub_folder, script)
    # cross_match在控制台输出序列映射到引物序列的映射关系
    # 去除引物py文件从控制台读取映射关系
    with open(cm_script, 'w') as output:
        for (ident, files) in split_file.items():
            for file in files:
                cmd = 'cross_match {} {}'.format(file, map_file)
                cmd += ' -gap1_only -minmatch 6 -minscore 10 -gap_init -3'
                cmd += ' | python3 {} {}'.format(script, file)
                cmd += '\n'
                output.write(cmd)
                cm_file = '{}.{}'.format(file, suffix)
                cm_files[ident].append(cm_file)
    return cm_script, cm_files


def main(work_dir, primer_file, idents, is_gzip, maxsize=200000):
    logging.info('STEP 2: REMOVING PRIMER SEQUENCES')
    tempdir = mkdtemp(dir=path.join(work_dir, 'temp'), prefix='primer_')
    script = 'CleanPrimer.py'
    # 将大文件分割为maxsize条序列的小文件，多进程同时进行任务
    suffix = 'fasta.gz' if is_gzip else 'fasta'
    split_file = read_fasta(work_dir, tempdir, idents, maxsize, suffix)
    script, result = make_script(split_file, primer_file, tempdir, script, suffix='clean')
    submit_array(script, 'primer_run', tempdir)
    # 合并清除序列后的分割文件
    files_list = merge_file(work_dir, result, suffix='clean.fasta')
    compress_multi(files_list.values())
    if Config()['auto_del']:
        rmtree(tempdir)
