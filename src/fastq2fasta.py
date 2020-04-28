#!/usr/bin/env python3
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from os import path, mkdir, remove, system
from random import randint
from subprocess import check_output

from src.singleton_config import Config


def fastq2fasta(fastq_file, fasta_file):
    change_cmd = 'seqtk seq -a {} > {}'.format(fastq_file, fasta_file)
    check_output(change_cmd, shell=True)
    logging.info('{} created'.format(path.basename(fasta_file)))


def merge_paired(fw_file, rc_file, fastq_file):
    merge_cmd = 'seqtk mergepe {} {} > {}'.format(fw_file, rc_file, fastq_file)
    check_output(merge_cmd, shell=True)
    logging.info('{} created'.format(path.basename(fastq_file)))


def random_select(fw_file, rc_file, count, filename, maxsize):
    # 随机筛选出maxsize条read，临时文件sample
    fw_sample = '{}.sample'.format(fw_file)
    rc_sample = '{}.sample'.format(rc_file)
    percent = maxsize / count
    # 使用相同的种子，使FW中随机筛选的read一一对应RC的read
    rand_seed = 100 + randint(100, 800)
    fw_seqtk = 'seqtk sample -s{} {} {} > {}'.format(rand_seed, fw_file, percent, fw_sample)
    rc_seqtk = 'seqtk sample -s{} {} {} > {}'.format(rand_seed, rc_file, percent, rc_sample)
    logging.info('There exists {} reads in the given fasta file. '
                 'Thus, approx. {:.2f} percent of reads will '
                 'be used in the remainder of assembly pipeline'.format(filename, percent))
    check_output(fw_seqtk, shell=True)
    check_output(rc_seqtk, shell=True)
    return fw_sample, rc_sample


def check_count(fw_file, rc_file, is_gzip):
    # fastq格式四行对应一条序列，每四行输出一行1，统计行数计算文件中的read数
    awk = "awk 'NR%4 == 1' {} | wc -l"
    if is_gzip:
        awk = "zcat {} | awk 'NR%4 == 1'| wc -l"
    fw_count = int(check_output(awk.format(fw_file), shell=True))
    rc_count = int(check_output(awk.format(rc_file), shell=True))
    return fw_count, rc_count


def compress(filename):
    cmd = 'gzip -f {}'.format(filename)
    system(cmd)
    logging.info('Gzip  {}'.format(filename))


def foo(work_dir, ident, fastq_dir, is_gzip, maxsize):
    out_file = path.join(work_dir, ident)
    # 创建样本的工作目录,后续中间结果文件均存放到对应的样本工作目录
    if not path.exists(out_file):
        mkdir(out_file)
    fw_file = path.join(fastq_dir, '{}_FW.fq'.format(ident))
    rc_file = path.join(fastq_dir, '{}_RC.fq'.format(ident))
    merge_file = path.join(fastq_dir, '{}.fastq'.format(ident))
    fasta_file = path.join(out_file, '{}.fasta'.format(ident))
    if is_gzip:
        fw_file += '.gz'
        rc_file += '.gz'
    # 判定分箱的fastq文件是否正确
    # 合并FW和RC文件后将fastq格式转换成fasta格式
    filename = path.basename(fw_file)
    fw_count, rc_count = check_count(fw_file, rc_file, is_gzip)
    logging.info('SampleID => {}\tFW_reads => {}\tRV_reads => {}'.format(filename, fw_count, rc_count))
    # FW的read数应该与RC中的read数相等
    # 如果样本的read数大于设定值，进行read筛选
    if fw_count == 0 and rc_count == 0:
        logging.info('No reads in {} file: 0 and 0'.format(filename))
    elif fw_count == rc_count:
        is_select = False
        if fw_count > maxsize:
            fw_file, rc_file = random_select(fw_file, rc_file, fw_count, filename, maxsize)
            is_select = True
        merge_paired(fw_file, rc_file, merge_file)
        fastq2fasta(merge_file, fasta_file)
        remove(merge_file)
        if is_select:
            remove(fw_file)
            remove(rc_file)
        if is_gzip:
            compress(fasta_file)
        return fasta_file
    else:
        logging.error('FW and RC file has unequal entries: {}: {} '
                      'and {}: {}'.format(fw_file, fw_count, rc_file, rc_count))


def main(work_dir, idents, is_gzip, maxsize):
    # 读取分箱后的fastq文件，每个样本对应一份FW和RC文件
    fastq_dir = path.join(work_dir, 'PreprocessedFiles')
    logging.info('Converting fastq to fasta format')
    logging.info('Maximum number of reads allowed:{}'.format(maxsize))
    if path.exists(fastq_dir):
        result = list()
        threads = Config()['thread']
        pool = ThreadPoolExecutor(threads)
        for ident in idents:
            result.append(pool.submit(foo, work_dir, ident, fastq_dir, is_gzip, maxsize))
        for future in as_completed(result):
            fasta_file = future.result()
            logging.info(fasta_file + 'had created.')
    else:
        logging.error('Error: Directory doesnt exists')
    logging.info("Format Conversion done successfully")
