#!/usr/bin/env python3

import sys
import fileinput
from re import search


def get_primer_seq():
    with fileinput.input(files='-') as stdin:
        seqs = set()
        primer = dict()
        switch = False
        # 从控制台读取序列映射到引物序列的映射关系
        for line in stdin:
            if 'Num. pairs' in line:
                switch = True
            if switch and line.startswith(' '):
                arr = [int(i) if i.isdigit() else i for i in line.split()]
                seq_name = arr[4]
                start_pos = arr[5]
                end_pos = arr[6]
                left_bp = int(search('(\d+)', arr[7]).group(0))
                # 引物映射的位置应该在序列的两端
                if start_pos <= 5:
                    primer[seq_name] = [0, end_pos]
                elif left_bp <= 5:
                    primer[seq_name] = [start_pos - 1, end_pos + left_bp]
                else:
                    seq_name = seq_name.split('#')[0]
                    seqs.add(seq_name)
            if search('matching entries', line):
                break
    return seqs, primer


def clean_primer(fasta_file, seqs, primer, out_file):
    with open(fasta_file) as file_in, open(out_file, 'w') as out:
        line = file_in.readline()
        while line:
            if line.startswith('>'):
                seq_name = line.lstrip('>').split('#')[0]

                name_fw = line.lstrip('>').strip()
                seq_fw = file_in.readline().strip()

                line = file_in.readline()
                name_rc = line.lstrip('>').strip()
                seq_rc = file_in.readline().strip()
                # 读取的序列必须为同一条
                if seq_name not in name_rc:
                    print("Second sequence doesn't match name")
                    exit(1)
                # N为序列中不确定的碱基，可表示为任意ACTG
                # 连续的N越长或N的数量越多，该序列质量越低
                if search('NNNNN+', seq_fw) \
                        or search('NNNN+', seq_rc) \
                        or seq_fw.count('N') + seq_rc.count('N') > 6 \
                        or seq_name in seqs:
                    line = file_in.readline()
                    continue
                # 移除引物序列对应的部分
                if name_fw in primer:
                    pos = primer[name_fw]
                    primer_seq = seq_fw[pos[0]:pos[1]]
                    seq_fw = seq_fw.replace(primer_seq, '')
                if name_rc in primer:
                    pos = primer[name_rc]
                    primer_seq = seq_rc[pos[0]:pos[1]]
                    seq_rc = seq_rc.replace(primer_seq, '')
                out.write('>{}\n'.format(name_fw))
                out.write(seq_fw + '\n')
                out.write('>{}\n'.format(name_rc))
                out.write(seq_rc + '\n')
            line = file_in.readline()


if __name__ == '__main__':
    fasta_file = sys.argv[1]
    seqs, primer = get_primer_seq()
    out_file = '{}.{}'.format(fasta_file, 'clean')
    clean_primer(fasta_file, seqs, primer, out_file)
