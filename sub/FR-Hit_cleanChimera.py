#!/usr/bin/env python3
import sys
from collections import defaultdict


def read_contigs(cd_file):
    length = defaultdict(int)
    contigs = defaultdict(str)
    with open(cd_file) as cd_in:
        for line in cd_in:
            if line.startswith('>'):
                name = line.lstrip('>').strip().strip('\n')
            else:
                contigs[name] += line.strip('\n')
                length[name] = len(contigs[name])
    return contigs, length


def create_mapped_list(seq_length):
    genomes = dict()
    for (name, length) in seq_length.items():
        genomes[name] = [0] * length
    return genomes


def get_mapping(fr_file, length):
    with open(fr_file) as fr_in:
        genomes_a = create_mapped_list(length)
        genomes_b = create_mapped_list(length)
        for line in fr_in:
            array = line.split()
            read_length = int(array[1].strip('nt'))
            alignment_length = int(array[3])
            identity = float(array[7].strip('%'))
            mapped_seq = array[8]
            mapped_start = int(array[9])
            mapped_end = int(array[10])

            if identity < 95:
                continue
            if not (alignment_length / read_length >= 0.95
                    or mapped_start < 5
                    or mapped_end > length[mapped_seq] - 5):
                continue
            for i in range(mapped_start, mapped_end + 1):
                genomes_a[mapped_seq][i - 1] += 1
                if mapped_start + 10 < i < mapped_end - 10:
                    genomes_b[mapped_seq][i - 1] += 1
                elif i <= 12 and i < mapped_end - 10:
                    genomes_b[mapped_seq][i - 1] += 1
                elif i > mapped_start + 10 and i > length[mapped_seq] - 15:
                    genomes_b[mapped_seq][i - 1] += 1
        return genomes_a, genomes_b


def clean_chimera(seqs, length, genomes_a, genomes_b, out_file):
    with open(out_file, 'w') as out:
        for (seq_name, coverage) in genomes_a.items():
            coverage_b = genomes_b[seq_name]
            coverage_b = [1 if i == 0 else i for i in coverage_b]
            for i in range(length[seq_name]):
                ratio = coverage[i] / coverage_b[i]
                if coverage_b[i] < 20 and ratio > 10 or ratio > 100:
                    seq = seqs[seq_name]
                    seqs[seq_name] = seq[:i] + 'N' + seq[i:]
            sub_seqs = seqs[seq_name].split('N')
            for i in range(len(sub_seqs)):
                if len(sub_seqs[i]) > 60:
                    sub_name = '>{}_{}\n'.format(seq_name, i)
                    sub_seq = sub_seqs[i] + '\n'
                    out.write(sub_name)
                    out.write(sub_seq)


if __name__ == '__main__':
    fr_file = sys.argv[1]
    seq_file = sys.argv[2]
    out_file = sys.argv[3]
    seqs, length = read_contigs(seq_file)
    genomes_a, genomes_b = get_mapping(fr_file, length)
    clean_chimera(seqs, length, genomes_a, genomes_b, out_file)
