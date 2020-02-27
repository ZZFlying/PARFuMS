#!/usr/bin/env python3


def read_diff():
    pl_put = '/home/ubuntu/wz4_fos.Missing1stPass.pl.fasta'
    py_put = '/home/ubuntu/wz4_fos.Missing1stPass.py.fasta'
    with open(pl_put) as pl_in, open(py_put) as py_in:
        py_seq = set()
        dif_seq = set()
        for line in py_in:
            if line.startswith('>'):
                py_seq.add(line)

        for line in pl_in:
            if line.startswith('>') and line not in py_seq:
                line = line.lstrip('>').strip('\n')
                dif_seq.add(line)
        return dif_seq


def read_sequence(dif_seq):
    fr_file = '/home/ubuntu/wz4_fos.Map2.txt'
    out_file = '/home/ubuntu/dir_map.txt'
    with open(fr_file) as fr_in, open(out_file, 'w') as out:
        for line in fr_in:
            seq_name = line.split()[0]
            if seq_name in dif_seq:
                out.write(line)


if __name__ == '__main__':
    dif_seq = read_diff()
    read_sequence(dif_seq)
