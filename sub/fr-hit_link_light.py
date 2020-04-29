#!/usr/bin/env python3
import sys
from os import getcwd, path, mkdir, system
from collections import defaultdict
from tempfile import TemporaryDirectory


def read_contigs(file):
    # 读取需要处理的序列
    contigs = defaultdict(str)
    length = defaultdict(int)
    with open(file) as cd_in:
        for line in cd_in:
            if line.startswith('>'):
                name = line.lstrip('>')
                name = name.strip()
            else:
                contigs[name] += line.strip()
                length[name] = len(contigs[name])
    return contigs, length


def process_array(array, length, links):
    map_seq = set()
    ref_seq = set()
    for mapping in array:
        map_name = mapping[8]
        ref_name = mapping[0]
        map_seq.add(map_name)
        ref_seq.add(ref_name)
    # 一条read的两条链各唯一映射一条长序列
    if len(array) != 2 or len(map_seq) != 2 or len(ref_seq) != 2:
        return

    seq_fw, seq_rc = array[0], array[1]
    seq_fw = [int(x) if x.isdigit() else x for x in seq_fw]
    seq_rc = [int(x) if x.isdigit() else x for x in seq_rc]
    map_fw, map_rc = seq_fw[8], seq_rc[8]
    # 生成两条长序列连接的方式
    # 长序列应该首尾相接
    if seq_fw[6] == '+':
        link = 'stop_1' if seq_fw[-1] > length[map_fw] - 100 else 'weird1'
    else:
        link = 'start_1' if seq_fw[-2] < 100 else 'weird_1'

    if seq_rc[6] == '+':
        link += '!stop_2' if seq_rc[-1] > length[map_rc] - 100 else '!weird2'
    else:
        link += '!start_2' if seq_rc[-2] < 100 else '!weird_2'
    # 统计两条长序列不同连接方式的次数
    if 'weird' not in link:
        fw_rc = '{}!{}'.format(map_fw, map_rc)
        rc_fw = '{}!{}'.format(map_rc, map_fw)
        if link == 'start_1!start_2' or link == 'stop_1!stop_2':
            if links[rc_fw][link]:
                links[rc_fw][link] += 1
            else:
                links[fw_rc][link] += 1
        elif link == 'start_1!stop_2':
            if links[rc_fw]['stop_1!start_2']:
                links[rc_fw]['stop_1!start_2'] += 1
            else:
                links[fw_rc][link] += 1
        elif link == 'stop_1!start_2':
            if links[rc_fw]['start_1!stop_2']:
                links[rc_fw]['start_1!stop_2'] += 1
            else:
                links[fw_rc][link] += 1


def check_links(fr_file, length):
    # 根据映射关系获得两条长序列的链接关系
    with open(fr_file) as fr_in:
        temp = ''
        mapping = list()
        links = defaultdict(lambda: defaultdict(int))
        for line in fr_in:
            array = line.split()
            seq_name = array[0].split('#')[0]
            coverage = float(array[7].strip('%'))
            if coverage < 95:
                continue
            if seq_name != temp:
                if temp != '':
                    process_array(mapping, length, links)
                temp = seq_name
                mapping = list()
            mapping.append(array)
        process_array(mapping, length, links)
        return links


def update_link_threshold(links):
    threshold = 5
    thresholds = [5, 10, 20, 40, 60, 80, 100]
    for threshold in thresholds:
        link_count = 0
        for link_seq in links.values():
            for count in link_seq.values():
                if count > threshold:
                    link_count += 1
        if link_count < 40:
            break
    return threshold


def create_tempdir(tempdir):
    directory = TemporaryDirectory(prefix='linkFiles', dir=tempdir)
    return directory


def get_sequences(link_map, fasta_file, out_file):
    with open(link_map) as link_in:
        seq_name = set()
        used = defaultdict(bool)
        for line in link_in:
            array = line.split()
            name = array[0].split('#')[0]
            if used[name]:
                continue
            if array[6] == '+' and 'start' in array[8]:
                continue
            if array[6] == '-' and 'stop' in array[8]:
                continue
            if 'stop' in array[8] and int(array[-1]) < 90:
                continue
            if 'start' in array[8] and int(array[-2]) > 10:
                continue
            used[name] = True
            seq_name.add(name + '#0_0')
            seq_name.add(name + '#0_1')

    with open(fasta_file) as fasta_in, open(out_file, 'a') as out:
        output = False
        for line in fasta_in:
            if line.startswith('>'):
                output = False
                name = line.lstrip('>')
                name = name.strip()
                if name in seq_name:
                    out.write(line)
                    output = True
            elif output:
                out.write(line)


def create_mapped_list(seq_length):
    genomes = dict()
    for (name, length) in seq_length.items():
        genomes[name] = [0] * length
    return genomes


def clean_chimera(fr_file, seq_file):
    seqs, length = read_contigs(seq_file)
    genomes_a = create_mapped_list(length)
    genomes_b = create_mapped_list(length)
    with open(fr_file) as fr_in:
        for line in fr_in:
            array = line.split()
            array[1] = array[1].strip('nt')
            array[7] = array[7].split('.')[0]
            array = [int(x) if x.isdigit() else x for x in array]
            if array[7] < 95:
                continue
            if not (array[3] / array[7] >= 0.95 or array[-2] < 5 or array[-1] > length[array[8]] - 5):
                continue
            for i in range(array[-2], array[-1] + 1):
                genomes_a[array[8]][i - 1] += 1
                if array[-2] + 10 < i < array[-1] - 10:
                    genomes_b[array[8]][i - 1] += 1
                elif i <= 12 and i < array[-1] - 10:
                    genomes_b[array[8]][i - 1] += 1
                elif i > array[-2] + 10 and i > length[array[8]] - 15:
                    genomes_b[array[8]][i - 1] += 1

    clean_seq = dict()
    for (seq_name, coverage) in genomes_a.items():
        coverage_b = genomes_b[seq_name]
        coverage_b = [1 if i == 0 else i for i in coverage_b]
        for i in range(length[seq_name]):
            precent = coverage[i] / coverage_b[i]
            if coverage_b[i] < 20 and precent > 10 or precent > 100:
                seq = seqs[seq_name]
                seqs[seq_name] = seq[:i] + 'N' + seq[i + 1:]
        sub_seqs = seqs[seq_name].split('N')
        for i in range(len(sub_seqs)):
            if len(sub_seqs[i]) > 60:
                clean_seq['{}_{}'.format(seq_name, i)] = sub_seqs[i]
    return clean_seq


def main(fr_file, cd_file, tempdir, ident, out_file):
    # 读取要处理的长序列
    contigs, length = read_contigs(cd_file)
    # 获得长序列的不同连接关系次数
    links = check_links(fr_file, length)
    # 确定进行处理的阈值
    threshold = update_link_threshold(links)

    dirname = path.dirname(out_file)
    fasta_file = path.join(dirname, ident + '.noVector.fasta')

    directory = create_tempdir(tempdir)
    tempdir = directory.name
    to_map_file = path.join(tempdir, "{}-Link_ToMap.fna".format(ident))
    link_map_file = path.join(tempdir, "{}-Link_Map.txt".format(ident))
    link_out_file = path.join(tempdir, "{}-Link_OutCdHit.fna".format(ident))

    count = 1
    modify = dict()
    clean_seq = dict()
    skip = defaultdict(bool)
    # 每对长序列的每种连接方式
    for (combine_name, link) in links.items():
        for (pos, link_num) in link.items():
            if link_num >= threshold:
                # 获取两条长序列名和连接方式
                name = combine_name.split('!')
                pos = pos.split('!')

                # 根据连接方式截取长序列的100个碱基，连接序列
                contig = contigs[name[0]]
                sub1 = contig[0:100] if 'start' in pos[0] else contig[-100:]
                contig = contigs[name[1]]
                sub2 = contig[0:100] if 'start' in pos[1] else contig[-100:]

                modify['{}.{}'.format(name[0], count)] = contigs[name[0]]
                modify['{}.{}'.format(name[1], count)] = contigs[name[1]]
                skip[name[0]] = True
                skip[name[1]] = True
                count += 1
                # 原始序列文件映射到连接序列，获取映射的序列
                # 追加到连接序列中进行组装
                # 输出去除嵌合体的序列
                with open(to_map_file, 'w') as output:
                    output.write('>{}_{}\n'.format(pos[0], name[0]))
                    output.write(sub1 + '\n')
                    output.write('>{}_{}\n'.format(pos[1], name[1]))
                    output.write(sub2 + '\n')
                system('fr-hit -d {} -a {} -o {} -g 1 -q 50 -c 95'.format(to_map_file, fasta_file, link_map_file))
                get_sequences(link_map_file, fasta_file, to_map_file)
                system('cd-hit-est -i {} -o {} -G 0 -aS 0.99 -g 1 -r 1 -c 0.9'.format(to_map_file, link_out_file))
                system("phrap -minmatch 10 -maxmatch 30 -bandwidth 0 -minscore 15 {}".format(link_out_file))
                system("fr-hit -d {}.contigs -o {} -a {} -m 30".format(link_out_file, link_map_file, fasta_file))
                # clean_seq.update(clean_chimera(link_map_file, link_out_file + '.contigs'))
                no_chimera = clean_chimera(link_map_file, link_out_file + '.contigs')
                for (k, v) in no_chimera.items():
                    clean_seq['{}.{}'.format(k, count)] = v
                    count += 1

    with open(out_file, 'w') as out:
        for (k, v) in contigs.items():
            if not skip[k]:
                out.write('>{}\n{}\n'.format(k, v))
        for (k, v) in modify.items():
            if len(v) > 30:
                out.write('>{}\n{}\n'.format(k, v))
        for (k, v) in clean_seq.items():
            if len(v) > 30:
                out.write('>{}\n{}\n'.format(k, v))

if __name__ == '__main__':
    [fr_file, cd_file, tempdir, ident, out_file] = sys.argv[1:6]
    main(fr_file, cd_file, tempdir, ident, out_file)
