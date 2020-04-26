#!/usr/bin/env python3
import sys
import logging
import traceback
from os import path, mkdir

import argparse

from src.singleton_config import Config
from src.load_barcode import load_barcode

from src.preprocessing import main as preprocessing
from src.fastq2fasta import main as fastq2fasta
from src.remove_primer import main as remove_primer
from src.remove_vector import main as remove_vector
from src.velvet_assembly import main as velvet_assembly
from src.phrap_assembly import main as phrap_assembly


def print_info():
    print('''
    # 双端序列文件，fastq.gz格式
    FW_file: /root/PARFuMS/datasets/FW.fastq.gz
    RC_file: /root/PARFuMS/datasets/RC.fastq.gz
    # Barcode文件
    # Barcode   Name
    # GGGGGGGG  样本名
    BC_file: /root/PARFuMS/examples/barcode.txt
    # 输出文件夹，按样本名区分
    work_dir: /root/PARFuMS/output
    # 引物和载体序列文件，fasta格式
    primer_file: /root/PARFuMS/datasets/primer.fasta
    vector_file: /root/PARFuMS/datasets/vector.fasta
    # 是否输出无法确定样本的read
    mismatch: true
    # 用于组装的最大read数
    maxsize: 8000000
    # 运行时的最大进程数
    thread: 8
    # 自动删除运行过程中产生的临时文件
    autodel: true
    # 日志记录等级
    # info = 0; debug = 1
    logger_level: 0
    ''')


def info():
    parse = argparse.ArgumentParser(description='Auto assembly NGS data.')
    parse.add_argument('config', help='配置文件路径', )
    parse.add_argument('step', help='运行指定步骤，为0全运行。'
                                    '多个步骤如第2步到第5步使用2:5表示。')
    _argv = parse.parse_args()
    return _argv

def init_logger(work_dir, logger_level):
    levels = [logging.INFO, logging.DEBUG]
    # 初始化日志,日志文件存放在工作目录下
    if logger_level > 1 or logger_level < 0:
        logger_level = 1
    logging.basicConfig(
        filename=path.join(work_dir, 'parfums.log'),
        level=levels[logger_level],
        format='%(asctime)s %(levelname)s %(module)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_pre_file(work_dir, idents):
    # 检测分箱的样本文件是否存在
    # 返回存在的样本标识
    fastq_dir = path.join(work_dir, 'PreprocessedFiles')
    for ident in idents:
        fw_file = path.join(fastq_dir, '{}_FW.fq'.format(ident))
        rc_file = path.join(fastq_dir, '{}_RC.fq'.format(ident))
        if path.exists(fw_file) and path.exists(rc_file):
            continue
        idents.remove(ident)
    return idents


def check_output_files(workdir, idents, step, suffix):
    # 检测主要步骤的结果文件
    # 如果文件不存在或文件大小为0,认为该步骤出现错误,移除对应的样本标识
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
    logging.debug('Check passed list:' + idents.__str__())
    return idents


def make_output_dir(work_dir):
    # 生成临时文件夹,用于存放过程中临时生成的样品分割文件
    preprocessed_files = path.join(work_dir, 'PreprocessedFiles')
    temp = path.join(work_dir, 'temp')
    if not path.exists(work_dir):
        mkdir(work_dir)
        mkdir(preprocessed_files)
        mkdir(temp)
    if not path.exists(preprocessed_files):
        mkdir(preprocessed_files)
    if not path.exists(temp):
        mkdir(temp)


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
        step_list.append(int(step))
    return step_list


if __name__ == '__main__':
    argv = info()
    step = argv.step
    config_file = argv.config
    # 生成单例的Config类, 方便其他模块直接访问参数
    config = Config(config_file)
    step_list = parse_step(step)

    work_dir = config['work_dir']
    primer_file = config['primer_file']
    vector_file = config['vector_file']
    fw_file = config['FW_file']
    rc_file = config['RC_file']
    bc_file = config['BC_file']
    maxsize = config['maxsize']
    mismatch = config['mismatch']
    logger_level = config['logger_level']

    # 生成工作目录，初始化日志，读取样本的barcode和标识
    make_output_dir(work_dir)
    init_logger(work_dir, logger_level)
    barcodes, idents, length = load_barcode(bc_file)

    logging.debug('Loaded barcodes: ' + barcodes.__str__())
    try:
        for step in step_list:
            if step == 1:
                # 样本分箱
                preprocessing(fw_file, rc_file, work_dir, barcodes, length, mismatch)
                idents = check_pre_file(work_dir, idents)
                # 样本read数筛选和格式转换
                fastq2fasta(work_dir, idents, maxsize)
                logging.info('STEP 1: PREPROCESSING COMPLETED')
            elif step == 2:
                # 去除引物
                idents = check_output_files(work_dir, idents, step, suffix='fasta')
                remove_primer(work_dir, primer_file, idents)
                logging.info('STEP 2: PRIMER SEQUENCES REMOVED')
            elif step == 3:
                # 去除载体
                idents = check_output_files(work_dir, idents, step, suffix='clean.fasta')
                remove_vector(work_dir, vector_file, idents)
                logging.info('STEP 3: VECTOR SEQUENCES REMOVED')
            elif step == 4:
                # velvet组装
                idents = check_output_files(work_dir, idents, step, suffix='noVector.fasta')
                velvet_assembly(work_dir, idents)
                logging.info('STEP 4: VELVET ASSEMBLY COMPLETED')
            elif step == 5:
                # phrap组装
                idents = check_output_files(work_dir, idents, step, suffix='ForPhrap1.fasta')
                phrap_assembly(work_dir, idents)
                logging.info('STEP 5: PHRAP ASSEMBLY COMPLETED')
            else:
                print('Step Count doesn\'t make sense in current program')
                logging.info('Step Count doesn\'t make sense in current program')
                exit(1)
    except KeyboardInterrupt:
        # 运行时强行退出
        logging.info('INTERRUPT THE PARFUMS PIPELINE')
    except Exception as e:
        # 运行中出现系统错误
        traceback.print_exc()
        logging.info('EXIT WITH ERROR')
        logging.info(traceback.format_exc())
    logging.info('END OF PARFUMS PIPELINE')

