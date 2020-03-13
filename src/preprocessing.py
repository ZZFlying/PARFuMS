#!/usr/bin/env python3
import logging
from os import path, mkdir


def create_output_files(work_dir, barcodes, mismatch):
    file_handles = dict()
    read_seqs = dict()
    tempdir = path.join(work_dir, 'PreprocessedFiles')
    if not path.exists(tempdir):
        mkdir(tempdir)
    # 创建每个样本对应FW和RC的文件io
    for (barcode, ident) in barcodes.items():
        filename = path.join(tempdir, ident)
        fw = open(filename + '_FW.fq', 'w')
        rc = open(filename + '_RC.fq', 'w')
        file_handles[ident] = fw, rc
        read_seqs[ident] = 0
    if mismatch:
        filename = path.join(tempdir, 'Mismatch')
        fw = open(filename + '_FW.fq', 'w')
        rc = open(filename + '_RC.fq', 'w')
        file_handles['Mismatch'] = fw, rc
        read_seqs['Mismatch'] = 0
    return file_handles, read_seqs


def write_summary(work_dir, read_seqs, total):
    match = 0
    summary_file = path.join(work_dir, 'SummaryFile.txt')
    with open(summary_file, 'w') as output:
        logging.info('Printing Summary to SummaryFile.txt')
        for (ident, count) in read_seqs.items():
            match += count
            record = '{}\t{}\n'.format(ident, count)
            output.write(record)
        mismatch = total - match
        mismatch = 'Mismatches\t{}\nTotal\t{}\n'.format(mismatch, total)
        output.write(mismatch)


def close_file(file_handles):
    for (file, handle) in file_handles.items():
        fw, rc = handle
        fw.close()
        rc.close()


def main(fw_file, rc_file, work_dir, barcodes, length, mismatch):
    try:
        with open(fw_file) as fw, open(rc_file) as rc:
            logging.info('STEP 1: PREPROCESSING OF FASTA FILES STARTED')

            fw_filename = path.basename(fw_file)
            rc_filename = path.basename(rc_file)
            logging.info('Reading sequence files:{} and {}'.format(fw_filename, rc_filename))

            total = 0
            fwd_line = fw.readline()
            rev_line = rc.readline()
            out_file_handles, read_seqs = create_output_files(work_dir, barcodes, mismatch)
            while fwd_line and rev_line:
                if fwd_line.startswith('@') and rev_line.startswith('@'):
                    total += 1

                    fwd_id = '{}#0_0\n'.format(fwd_line.split()[0])
                    fwd_bar = fwd_line.strip()[-length:]
                    fwd_seq = fw.readline()
                    fwd_plus = fw.readline()
                    fwd_qual = fw.readline()

                    rev_id = '{}#0_1\n'.format(rev_line.split()[0])
                    rev_bar = rev_line.strip()[-length:]
                    rev_seq = rc.readline()
                    rev_plus = rc.readline()
                    rev_qual = rc.readline()
                    # 根据读取的barcode获得对应样本的文件io
                    # 输出读取的序列到对应的文件
                    if fwd_bar in barcodes and fwd_bar == rev_bar:
                        out_fw, out_rc = out_file_handles[barcodes[fwd_bar]]
                        out_fw.writelines([fwd_id, fwd_seq, fwd_plus, fwd_qual])
                        out_rc.writelines([rev_id, rev_seq, rev_plus, rev_qual])
                        read_seqs[barcodes[fwd_bar]] += 1
                    elif mismatch:
                        # 输出无法匹配任何样本的序列
                        out_fw, out_rc = out_file_handles['Mismatch']
                        out_fw.writelines([fwd_id, fwd_seq, fwd_plus, fwd_qual])
                        out_rc.writelines([rev_id, rev_seq, rev_plus, rev_qual])
                        read_seqs['Mismatch'] += 1

                fwd_line = fw.readline()
                rev_line = rc.readline()
            # 安全关闭文件io，保障所有序列的正常写入文件
            close_file(out_file_handles)
            write_summary(work_dir, read_seqs, total)
    except FileNotFoundError:
        logging.error('Failed to open file')
