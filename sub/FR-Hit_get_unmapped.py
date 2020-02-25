#!/usr/bin/env python3
import sys


def process_array_old(arrays, used):
    seq_name = arrays[0][0].split('#')[0]
    read_length = arrays[0][1]
    count = 0
    # 原始seq可能映射到多个contig序列，合并各个映射区间
    seq = [0] * (read_length + 1)
    for arr in arrays:
        [start, end] = sorted(arr[4:6])
        for i in range(start, end + 1):
            seq[i] += 1
    # 计算合并区间的长度，超过80%被映射，记该seq被使用
    for i in range(1, read_length + 1):
        if seq[i]:
            count += 1
    if count > read_length * 0.8:
        used.add(seq_name)


def process_array(arrays, used):
    count = 0
    # 原始seq可能映射到多个contig序列，合并各个映射区间
    interval = list()
    for arr in arrays:
        interval.append(sorted([arr[4], arr[5]]))
    interval.sort()
    res = list()
    for i in interval:
        if not res or res[-1][1] < i[0]:
            res.append(i)
        else:
            res[-1][1] = max(res[-1][1], i[1])
    # 计算映射的总长度
    for i in res:
        count += i[1] - i[0] + 1
    seq_name = arrays[0][0].split('#')[0]
    read_length = arrays[0][1]
    if count > read_length * 0.8:
        used.add(seq_name)


def get_seq_name(fr_file):
    with open(fr_file) as fr_in:
        temp = ''
        seq_list = list()
        used = set()
        for line in fr_in:
            array = line.split()
            array[1] = array[1].strip('nt')
            array = [int(i) if i.isdigit() else i for i in array]
            seq_name = array[0]
            if seq_name != temp:
                if temp:
                    process_array(seq_list, used)
                temp = seq_name
                seq_list = list()
            seq_list.append(array)
        process_array(seq_list, used)
        seq_name = set()
        for n in used:
            seq_fw = n + '#0_0'
            seq_rc = n + '#0_1'
            seq_name.add(seq_fw)
            seq_name.add(seq_rc)
        return seq_name


def get_unmapped_seq(seq_file, seq_name, out_file):
    with open(seq_file) as seq_in, open(out_file, 'w') as out:
        for line in seq_in:
            if line.startswith('>'):
                name = line.lstrip('>').strip()
                if name not in seq_name:
                    out.write(line)
                    writeable = True
                else:
                    writeable = False
            elif writeable:
                out.write(line)


if __name__ == '__main__':
    [fr_file, seq_file, out_file] = sys.argv[1:4]
    seq_name = get_seq_name(fr_file)
    get_unmapped_seq(seq_file, seq_name, out_file)

