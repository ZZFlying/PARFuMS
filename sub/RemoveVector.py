#!/usr/bin/env python3

import sys
import fileinput
from re import search


def get_vector_seq():
    # 从控制台读取序列映射到引物序列的映射关系
    with fileinput.input(files='-') as stdin:
        switch = False
        vector_seq = set()
        multi_ref = set()
        for line in stdin:
            if 'Num. pairs' in line:
                switch = True
            if switch and line.startswith(' '):
                arr = [int(i) if i.isdigit() else i for i in line.split()]
                seq_name = arr[4]
                start_pos = arr[5]
                end_pos = arr[6]
                left_bp = int(arr[7][1:-1])
                is_vector = False
                # 多次映射到载体序列
                if seq_name in multi_ref:
                    is_vector = True
                # 映射在序列两侧，或覆盖度大于85%，丢弃该序列
                if start_pos <= 5 or left_bp <= 5:
                    is_vector = True
                if abs((end_pos - start_pos) / (end_pos + left_bp)) > 0.85:
                    is_vector = True
                multi_ref.add(seq_name)
                if is_vector:
                    paired_name = seq_name.split('#')[0]
                    vector_seq.add(paired_name)
            if search('matching entries', line):
                break
    return vector_seq


def clean_vector(fasta_file, vector, out_file):
    with open(fasta_file) as file_in, open(out_file, 'w') as out:
        line = file_in.readline()
        while line:
            if line.startswith('>'):
                seq_name = line.lstrip('>').split('#')[0]

                name_fw = line
                seq_fw = file_in.readline()

                name_rc = file_in.readline()
                seq_rc = file_in.readline()
                if seq_name not in name_rc:
                    print("Second sequence doesn't match name")
                    exit(1)

                if seq_name not in vector:
                    out.write(name_fw)
                    out.write(seq_fw)
                    out.write(name_rc)
                    out.write(seq_rc)
            line = file_in.readline()


if __name__ == '__main__':
    fasta_file = sys.argv[1]
    out_file = '{}.{}'.format(fasta_file, 'noVector')
    vector = get_vector_seq()
    clean_vector(fasta_file, vector, out_file)
