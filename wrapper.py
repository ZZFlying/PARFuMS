#!/usr/bin/env python3

import logging
import traceback
import argparse
from os import path, mkdir

from util.singleton_config import Config
from util.load_barcode import load_barcode

from src.preprocessing import main as preprocessing
from src.fastq2fasta import main as fastq2fasta
from src.remove_primer import main as remove_primer
from src.remove_vector import main as remove_vector
from src.velvet_assembly import main as velvet_assembly
from src.phrap_assembly import main as phrap_assembly


def info():
    parse = argparse.ArgumentParser(description='Auto assembly NGS data.')
    parse.add_argument('config', help='配置文件路径', )
    parse.add_argument('step', help='运行指定步骤，为0全运行。'
                                    '多个步骤如第2步到第5步使用2:5表示。')
    _argv = parse.parse_args()
    return _argv


def init_logger(_work_dir, _logger_level):
    levels = [logging.INFO, logging.DEBUG]
    # 初始化日志,日志文件存放在工作目录下
    if _logger_level > 1 or _logger_level < 0:
        _logger_level = 1
    logging.basicConfig(
        filename=path.join(_work_dir, 'parfums.log'),
        level=levels[logger_level],
        format='%(asctime)s %(levelname)s %(module)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def check_pre_file(_work_dir, _idents):
    # 检测分箱的样本文件是否存在
    # 返回存在的样本标识
    fastq_dir = path.join(_work_dir, 'PreprocessedFiles')
    for ident in _idents:
        fw = path.join(fastq_dir, '{}_FW.fq'.format(ident))
        rc = path.join(fastq_dir, '{}_RC.fq'.format(ident))
        if path.exists(fw) and path.exists(rc):
            continue
        if path.exists(fw + '.gz') and path.exists(rc + '.gz'):
            continue
        _idents.remove(ident)
    return _idents


def check_output_files(_word_dir, _idents, _step, _suffix):
    # 检测主要步骤的结果文件
    # 如果文件不存在或文件大小为0,认为该步骤出现错误,移除对应的样本标识
    not_pass_list = list()
    for ident in _idents:
        out_file = path.join(_word_dir, ident, '{}.{}'.format(ident, _suffix))
        if path.exists(out_file) and path.getsize(out_file):
            continue
        logging.info('{} is empty or not found. Halting its processing here'.format(out_file))
        not_pass_list.append(ident)
        _idents.remove(ident)
    if not_pass_list:
        not_pass_seq = ' '.join([i for i in not_pass_list])
        logging.info('IDs not included at Step{}: {}'.format(_step, not_pass_seq))
    logging.debug('Check passed list:' + _idents.__str__())
    return _idents


def make_output_dir(_work_dir):
    # 生成临时文件夹,用于存放过程中临时生成的样品分割文件
    preprocessed_files = path.join(_work_dir, 'PreprocessedFiles')
    temp = path.join(_work_dir, 'temp')
    if not path.exists(_work_dir):
        mkdir(_work_dir)
        mkdir(preprocessed_files)
        mkdir(temp)
    if not path.exists(preprocessed_files):
        mkdir(preprocessed_files)
    if not path.exists(temp):
        mkdir(temp)


def parse_step(_step):
    _step = _step.strip()
    _step_list = list()
    if ':' in _step:
        _step = _step.split(':')
        start = int(step[0])
        stop = int(step[1])
        _step_list = range(start, stop + 1)
    elif int(step) == 0:
        _step_list = range(1, 6)
    else:
        _step_list.append(int(_step))
    return _step_list


if __name__ == '__main__':
    argv = info()
    step = argv.step
    config_file = argv.config
    # 生成单例的Config类, 方便其他模块直接访问参数
    config = Config(config_file)
    step_list = parse_step(step)

    is_gzip = config['gzip']
    fw_file = config['FW_file']
    rc_file = config['RC_file']
    bc_file = config['BC_file']
    maxsize = config['maxsize']
    work_dir = config['work_dir']
    mismatch = config['mismatch']
    primer_file = config['primer_file']
    vector_file = config['vector_file']
    logger_level = config['logger_level']
    # 生成工作目录，初始化日志，读取样本的barcode和标识
    make_output_dir(work_dir)
    init_logger(work_dir, logger_level)
    barcodes, idents = load_barcode(bc_file)

    logging.debug('Loaded barcodes: ' + barcodes.__str__())
    try:
        for step in step_list:
            if step == 1:
                # 样本分箱
                preprocessing(fw_file, rc_file, work_dir, barcodes, mismatch)
                idents = check_pre_file(work_dir, idents)
                # 样本read数筛选和格式转换
                fastq2fasta(work_dir, idents, is_gzip, maxsize)
                logging.info('STEP 1: PREPROCESSING COMPLETED')
            elif step == 2:
                # 去除引物
                suffix = 'fasta.gz' if is_gzip else 'fasta'
                idents = check_output_files(work_dir, idents, step, suffix)
                remove_primer(work_dir, primer_file, idents, is_gzip)
                logging.info('STEP 2: PRIMER SEQUENCES REMOVED')
            elif step == 3:
                # 去除载体
                idents = check_output_files(work_dir, idents, step, 'clean.fasta')
                remove_vector(work_dir, vector_file, idents)
                logging.info('STEP 3: VECTOR SEQUENCES REMOVED')
            elif step == 4:
                # velvet组装
                idents = check_output_files(work_dir, idents, step, 'noVector.fasta')
                velvet_assembly(work_dir, idents)
                logging.info('STEP 4: VELVET ASSEMBLY COMPLETED')
            elif step == 5:
                # phrap组装
                idents = check_output_files(work_dir, idents, step, 'ForPhrap1.fasta')
                phrap_assembly(work_dir, idents)
                logging.info('STEP 5: PHRAP ASSEMBLY COMPLETED')
            else:
                print("Step Count doesn't make sense in current program")
                logging.info("Step Count doesn't make sense in current program")
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

