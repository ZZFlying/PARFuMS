#!/usr/bin/env python3
import logging
from os import path, mkdir
from random import randint
from subprocess import check_output, check_call, CalledProcessError


def merge_fastq(fw_file, rc_file, merge_file):
    # 合并FW和RC的fastq文件，FW的read后紧跟RC对应的read
    # @seq#0_0
    # ...
    # @seq#0_1
    # ...
    with open(fw_file) as fw, open(rc_file) as rc, open(merge_file, 'w') as output:
        fwd_line = fw.readline()
        rev_line = rc.readline()
        while fwd_line and rev_line:
            if fwd_line.startswith('@') and rev_line.startswith('@'):
                fwd_id = fwd_line.split('#')[0]
                rev_id = rev_line.split('#')[0]
                if fwd_id == rev_id:
                    fwd_seq = fw.readline()
                    fwd_plus = fw.readline()
                    fwd_qual = fw.readline()

                    rev_seq = rc.readline()
                    rev_plus = rc.readline()
                    rev_qual = rc.readline()
                    output.writelines([fwd_line, fwd_seq, fwd_plus, fwd_qual])
                    output.writelines([rev_line, rev_seq, rev_plus, rev_qual])
            fwd_line = fw.readline()
            rev_line = rc.readline()


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
    check_call(fw_seqtk, shell=True)
    check_call(rc_seqtk, shell=True)
    return fw_sample, rc_sample


def check_read(fw_file, rc_file, merge_file, maxsize):
    # fastq格式四行对应一条序列，每四行输出一行1，统计行数计算文件中的read数
    awk = "awk 'NR%4 == 1' {} | wc -l"
    fw_count = check_output(awk.format(fw_file), shell=True)
    fw_count = int(fw_count)
    rc_count = check_output(awk.format(rc_file), shell=True)
    rc_count = int(rc_count)
    filename = path.basename(fw_file)
    logging.info('SampleID => {}\tFW_reads => {}\tRV_reads => {}'.format(filename, fw_count, rc_count))
    # FW的read数应该与RC中的read数相等
    # 如果样本的read数大于设定值，进行read筛选
    if fw_count == 0 and rc_count == 0:
        logging.info('No reads in {} file: 0 and 0'.format(filename))
    elif fw_count == rc_count:
        if fw_count > maxsize:
            fw_file, rc_file = random_select(fw_file, rc_file, fw_count, filename, maxsize)
        merge_fastq(fw_file, rc_file, merge_file)
        return True
    else:
        logging.error('FW and RV file has unequal entries: {}: {} '
                      'and {}: {}'.format(fw_file, fw_count, rc_file, rc_count))
    return False


def main(work_dir, idents, maxsize):
    # 读取分箱后的fastq文件，每个样本对应一份FW和RC文件
    fastq_dir = path.join(work_dir, 'PreprocessedFiles')
    logging.info('Converting fastq to fasta format')
    logging.info('Maximum number of reads allowed:{}'.format(maxsize))
    if path.exists(fastq_dir):
        for ident in idents:
            out_file = path.join(work_dir, ident)
            # 创建样本的工作目录,后续中间结果文件均存放到对应的样本工作目录
            if not path.exists(out_file):
                mkdir(out_file)
            fw_file = path.join(fastq_dir, '{}_FW.fq'.format(ident))
            rc_file = path.join(fastq_dir, '{}_RC.fq'.format(ident))
            merge_file = path.join(fastq_dir, '{}.fastq'.format(ident))
            fasta_file = path.join(out_file, '{}.fasta'.format(ident))
            # 判定分箱的fastq文件是否正确
            # 合并FW和RC文件后将fastq格式转换成fasta格式
            if check_read(fw_file, rc_file, merge_file, maxsize):
                fq2fa = 'seqtk seq -a {} > {}'.format(merge_file, fasta_file)
                try:
                    check_call(fq2fa, shell=True)
                    logging.info('{} created'.format(path.basename(fasta_file)))
                except CalledProcessError:
                    logging.error('Error in executing seq-tk program: {}'.format(fq2fa))
    else:
        logging.error('Error: Directory doesnt exists')
