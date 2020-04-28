#!/usr/bin/env python3
import logging
from os import path, mkdir, system
from concurrent.futures import ThreadPoolExecutor, as_completed

from src.fastqReader import fastqReader
from src.singleton_config import Config


def create_handle(file_handles, read_count, tempdir, barcode='Mismatch', ident='Mismatch'):
    filename = path.join(tempdir, ident + '_{}.fq')
    fw_file = filename.format('FW')
    rc_file = filename.format('RC')
    fw = open(fw_file, 'wt')
    rc = open(rc_file, 'wt')
    file_handles[barcode] = (fw, rc)
    read_count[barcode] = 0
    return [fw_file, rc_file]


def create_output_files(work_dir, barcodes, mismatch):
    files_list = list()
    read_count = dict()
    file_handles = dict()
    tempdir = path.join(work_dir, 'PreprocessedFiles')
    if not path.exists(tempdir):
        mkdir(tempdir)
    # 创建每个样本对应FW和RC的文件io
    for (barcode, ident) in barcodes.items():
        files_list.extend(create_handle(file_handles, read_count, tempdir, barcode, ident))
    if mismatch:
        files_list.extend(create_handle(file_handles, read_count, tempdir))
    return file_handles, read_count, files_list


def write_summary(work_dir, read_seqs, total):
    match = 0
    summary_file = path.join(work_dir, 'SummaryFile.txt')
    with open(summary_file, 'w') as output:
        logging.info('Printing Summary to SummaryFile.txt')
        for (ident, count) in read_seqs.items():
            match += count
            record = '{}    {}\n'.format(ident, count)
            output.write(record)
        mismatch = total - match
        mismatch = 'Mismatches  {}\nTotal   {}\n'.format(mismatch, total)
        output.write(mismatch)


def close_file(file_handles):
    for (file, handle) in file_handles.items():
        (fw, rc) = handle
        fw.close()
        rc.close()


def write_record(handles, barcode, fw_read, rc_read):
    (fw_handle, rc_handle) = handles[barcode]
    fw_handle.write(fw_read.parse(0))
    rc_handle.write(rc_read.parse(1))


def compress(filename):
    cmd = 'gzip -f {}'.format(filename)
    system(cmd)
    return filename


def compress_multi(files_list):
    threads = Config()['thread']
    pool = ThreadPoolExecutor(threads)
    result = list()
    for file in files_list:
        logging.info('Gzip  {}'.format(file))
        result.append(pool.submit(compress, file))
    for future in as_completed(result):
        gz_file = future.result()
        gz_file = gz_file + '.gz'
        logging.info(gz_file + ' had created.')


def main(fw_file, rc_file, work_dir, barcodes, is_gzip, mismatch):
    try:
        logging.info('STEP 1: PREPROCESSING OF FASTA FILES STARTED')
        with fastqReader(fw_file) as fw_in, fastqReader(rc_file) as rc_in:
            total = 0
            fw_filename = path.basename(fw_file)
            rc_filename = path.basename(rc_file)
            file_handles, read_count, files_list = create_output_files(work_dir, barcodes, mismatch)
            logging.info('Reading sequence files:{} and {}'.format(fw_filename, rc_filename))
            for (fw_read, rc_read) in zip(fw_in, rc_in):
                total += 1
                barcode = fw_read.barcode
                if barcode in barcodes and fw_read.barcode == rc_read.barcode:
                    write_record(file_handles, barcode, fw_read, rc_read)
                    read_count[barcode] += 1
                elif mismatch:
                    # 输出无法匹配任何样本的序列
                    write_record(file_handles, 'Mismatch', fw_read, rc_read)
                    read_count['Mismatch'] += 1
            # 安全关闭文件io，保障所有序列的正常写入文件
            close_file(file_handles)
            write_summary(work_dir, read_count, total)
            if is_gzip:
                compress_multi(files_list)
    except FileNotFoundError:
        logging.error('Failed to open file')
