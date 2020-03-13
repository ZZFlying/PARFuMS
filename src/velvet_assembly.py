#!/usr/bin/env python3
import sys
import logging
from os import path, mkdir

from sub.parfums_subs import create_tempdir, submit_array, read_fasta

step = 0


def read_inputSequences(work_dir, idents, suffix):
    # 返回suffix后缀文件的绝对路径
    input_sequences = dict()
    for ident in idents:
        filename = path.join(work_dir, ident, '{}.{}'.format(ident, suffix))
        if path.exists(filename):
            input_sequences[ident] = filename
        else:
            logging.error('{} not exist'.format(filename))
    logging.debug('read_inputSequences => ' + input_sequences.__str__())
    logging.info('Seq file stored')
    return input_sequences


def merge_velvet(work_dir, temp_dir, split_file, suffix):
    contigs = dict()
    for (ident, files) in split_file.items():
        out_file = path.join(work_dir, ident, '{}.{}'.format(ident, suffix))
        logging.info('Now merging {} velvet contigs to {}'.format(ident, out_file))
        with open(out_file, 'w') as out:
            for file in files:
                split_name = path.basename(file).split('.')[0]
                contig_file = path.join(temp_dir, 'assemble_{}'.format(split_name), 'contigs.fa')
                with open(contig_file) as contig_in:
                    out.writelines(contig_in.readlines())
            contigs[ident] = out_file
            logging.info('{} contigs merge completely'.format(out_file))
    return contigs


def get_velvet(split_file, temp_dir, params):
    script = list()
    (velveth_param, velvetg_param) = params
    for (ident, files) in split_file.items():
        for file in files:
            filename = path.basename(file).split('.')[0]
            assemble_dir = path.join(temp_dir, 'assemble_' + filename)
            velveth_cmd = ' '.join(['velveth', assemble_dir, velveth_param, file])
            velvetg_cmd = ' '.join(['velvetg', assemble_dir, velvetg_param])
            temp = '{} && {}\n'.format(velveth_cmd, velvetg_cmd)
            logging.debug('cmd => ' + temp)
            script.append(temp)
    return script


def get_cd_hit_est(contig_file, suffix):
    script = list()
    cdhit_file = dict()
    for (ident, file) in contig_file.items():
        dirname = path.dirname(file)
        out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
        cmd = 'cd-hit-est -i {} -o {} -g 1 -r 1 -c 0.9\n'.format(file, out_file)
        logging.debug('cmd => ' + cmd)
        script.append(cmd)
        cdhit_file[ident] = out_file
    return script, cdhit_file


def get_fr_hit(seq_file, cdhit_file, suffix):
    script = list()
    frhit_file = dict()
    for (ident, mapping_file) in seq_file.items():
        dirname = path.dirname(mapping_file)
        out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
        cmd = 'fr-hit -a {} -d {} -m 30 -o {}\n'.format(mapping_file, cdhit_file[ident], out_file)
        logging.debug('cmd => ' + cmd)
        script.append(cmd)
        frhit_file[ident] = out_file
    return script, frhit_file


def get_remove_chimera(cdhit_file, frhit_file, idents, suffix):
    script = list()
    no_chimera_file = dict()
    perl_script = sys.path[0]
    perl_script = path.join(perl_script, 'sub', 'FR-Hit_cleanChimera.py')
    for ident in idents:
        fr_file = frhit_file[ident]
        cd_file = cdhit_file[ident]
        dirname = path.dirname(cd_file)
        out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
        cmd = str.join(' ', ['python3', perl_script, fr_file, cd_file, out_file])
        cmd = cmd + '\n'
        logging.debug('cmd => ' + cmd)
        script.append(cmd)
        no_chimera_file[ident] = out_file
    return script, no_chimera_file


def get_unmapped_reads(seq_file, frhit_file, split_file, suffix):
    script = list()
    perl_script = sys.path[0]
    perl_script = path.join(perl_script, 'sub', 'FR-Hit_get_unmapped.py')
    # 获取没有被使用的read序列
    for (ident, file) in seq_file.items():
        fr_file = frhit_file[ident]
        dirname = path.dirname(file)
        out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
        cmd = str.join(' ', ['python3', perl_script, fr_file, file, out_file])
        cmd = cmd + '\n'
        logging.debug('cmd => ' + cmd)
        script.append(cmd)
        split_file[ident] = [out_file]
    return script, split_file


def get_combine_contig(idents, cdhit_file, contig_file, suffix):
    script = list()
    for ident in idents:
        dirname = path.dirname(cdhit_file[ident])
        out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
        cmd = str.join(' ', ['cat', cdhit_file[ident], contig_file[ident], '>', out_file])
        cmd = cmd + '\n'
        logging.debug('cmd => ' + cmd)
        contig_file[ident] = out_file
        script.append(cmd)
    return script


def make_script(run_type, temp_dir, suffix=None, idents=None, params=None,
                split=None, contig=None, seq=None, cdhit=None, frhit=None):
    global step
    step += 1
    job_filename = '{}_script_V{}.sh'.format(run_type, step)
    script_dir = path.join(temp_dir, 'script')
    if not path.exists(script_dir):
        mkdir(script_dir)
    job_file = path.join(script_dir, job_filename)
    script = list()
    out_file = dict()
    with open(job_file, 'w') as out:
        logging.info('Generating script file to run {}:{}'.format(run_type, job_filename))
        if run_type == 'velvet':
            logging.info('Writing Velvet script file')
            script = get_velvet(split, temp_dir, params)
        elif run_type == 'cd-hit-est':
            logging.info('Writing CD-HIT script file')
            script, out_file = get_cd_hit_est(contig, suffix)
        elif run_type == 'fr-hit':
            logging.info('Writing FR-HIT script file')
            script, out_file = get_fr_hit(seq, cdhit, suffix)
        elif run_type == 'remove-chimera':
            logging.info('Writing Remove-Chimera script file')
            script, out_file = get_remove_chimera(cdhit, frhit, idents, suffix)
        elif run_type == 'unmapped-reads':
            logging.info('Writing script file to get unmapped reads')
            script, out_file = get_unmapped_reads(seq, frhit, split, suffix)
        elif run_type == 'combine-contigs':
            logging.info('Writing script to combine contig files')
            script = get_combine_contig(idents, cdhit, contig, suffix)
        else:
            logging.error('Invalid runType:{}'.format(run_type))

        if len(script) == 0:
            logging.info('Nothing in the command list to submit on cluster!')
        else:
            out.writelines(script)
    return job_file, out_file


def round_1(work_dir, seq_file, idents):
    logging.info('Velvet Assembly Round-1 Started')
    directory = create_tempdir(work_dir, prefix='VelvetRun1_')
    temp_dir = directory.name

    velveth_param = '31 -shortPaired'
    velvetg_param = '-cov_cutoff 10 -ins_length 100 -min_contig_lgth 100'
    split_file = read_fasta(work_dir, temp_dir, idents, 10000, suffix='noVector.fasta')

    script, empty_file = make_script('velvet', temp_dir, split=split_file,
                                     params=(velveth_param, velvetg_param))
    submit_array(script, 'VelvetRun1', temp_dir)
    contig_file = merge_velvet(work_dir, temp_dir, split_file, 'velvet1_contigs.fasta')

    script, cdhit_file = make_script('cd-hit-est', temp_dir, contig=contig_file, suffix='cd-hit1.fasta')
    submit_array(script, 'CD-hit1', temp_dir)

    script, frhit_file = make_script('fr-hit', temp_dir, seq=seq_file, cdhit=cdhit_file, suffix='Map1.txt')
    submit_array(script, 'FR-hit1', temp_dir)

    script, no_chimera_file = make_script('remove-chimera', temp_dir, cdhit=cdhit_file, frhit=frhit_file,
                                          idents=idents, suffix='NoChimera.fasta')
    submit_array(script, 'Remove-Chimera1', temp_dir)

    script, cdhit_file = make_script('cd-hit-est', temp_dir, contig=no_chimera_file, suffix='cd-hit2.fasta')
    submit_array(script, 'CD-hit2', temp_dir)
    # 将原始的seq片段映射到Velvet组装后去重的contig片段
    script, frhit_file = make_script('fr-hit', temp_dir, seq=seq_file, cdhit=cdhit_file, suffix='Map2.txt')
    submit_array(script, 'FR-hit2', temp_dir)
    # 获取Fr-Hit中seq映射到contig的序列名
    # 读取原始seq片段，输出未映射的seq序列，存储到split_file，进行下一轮Velvet组装
    script, split_file = make_script('unmapped-reads', temp_dir,
                                     seq=seq_file, frhit=frhit_file, split=split_file, suffix='Missing1stPass.fasta')
    submit_array(script, 'UnmappedReads1', temp_dir)

    logging.info('Round-1 completed')
    return split_file, cdhit_file


def round_2(work_dir, seq_file, split_file, cdhit_file, idents):
    logging.info('Velvet Assembly Round-2 Started')
    directory = create_tempdir(work_dir, prefix='VelvetRun2_')
    temp_dir = directory.name

    velveth_param = '31 -shortPaired'
    velvetg_param = '-cov_cutoff 7 -ins_length 80 -min_contig_lgth 100'

    script, empty_file = make_script('velvet', temp_dir, split=split_file, params=(velveth_param, velvetg_param))
    submit_array(script, 'VelvetRun2', temp_dir)
    contig_file = merge_velvet(work_dir, temp_dir, split_file, 'velvet2_contigs.fasta')

    script, empty_file = make_script('combine-contigs', temp_dir,
                                     idents=idents, cdhit=cdhit_file, contig=contig_file, suffix='ForCD-hit.fasta')
    submit_array(script, 'CombineContigs1', temp_dir)

    script, cdhit_file = make_script('cd-hit-est', temp_dir, contig=contig_file, suffix='cd-hit3.fasta')
    submit_array(script, 'CD-hit3', temp_dir)

    script, frhit_file = make_script('fr-hit', temp_dir, seq=seq_file, cdhit=cdhit_file, suffix='Map3.txt')
    submit_array(script, 'FR-hit3', temp_dir)

    script, empty_file = make_script('unmapped-reads', temp_dir,
                                     seq=seq_file, frhit=frhit_file, split=split_file, suffix='Missing2ndPass.fasta')
    submit_array(script, 'UnmappedReads2', temp_dir)

    logging.info('Round-2 completed')
    return split_file, cdhit_file


def round_3(work_dir, seq_file, split_file, cdhit_file, idents):
    logging.info('Velvet Assembly Round-3 Started')
    directory = create_tempdir(work_dir, prefix='VelvetRun3_')
    temp_dir = directory.name

    velveth_param = '31 -shortPaired'
    velvetg_param = '-cov_cutoff 10 -ins_length 80 -min_contig_lgth 100'

    script, empty_file = make_script('velvet', temp_dir, split=split_file, params=(velveth_param, velvetg_param))
    submit_array(script, 'VelvetRun3', temp_dir)
    contig_file = merge_velvet(work_dir, temp_dir, split_file, suffix='velvet3_contigs.fasta')

    script, empty_file = make_script('combine-contigs', temp_dir,
                                     idents=idents, cdhit=cdhit_file, contig=contig_file, suffix='ForCD-hit2.fasta')
    submit_array(script, 'CombineContigs2', temp_dir)

    script, cdhit_file = make_script('cd-hit-est', temp_dir, contig=contig_file, suffix='cd-hit4.fasta')
    submit_array(script, 'CD-hit4', temp_dir)

    script, frhit_file = make_script('fr-hit', temp_dir, seq=seq_file, cdhit=cdhit_file, suffix='Map4.txt')
    submit_array(script, 'FR-hit4', temp_dir)

    script, empty_file = make_script('remove-chimera', temp_dir,
                                     cdhit=cdhit_file, frhit=frhit_file, idents=idents, suffix='ForPhrap1.fasta')
    submit_array(script, 'Remove-Chimera2', temp_dir)
    logging.info('Round-3 completed')


def main(work_dir, idents):
    logging.info('STEP 4: VELVET ASSEMBLY STARTED')

    seq_file = read_inputSequences(work_dir, idents, 'noVector.fasta')
    split_file, cdhit_file = round_1(work_dir, seq_file, idents)
    split_file, cdhit_file = round_2(work_dir, seq_file, split_file, cdhit_file, idents)
    round_3(work_dir, seq_file, split_file, cdhit_file, idents)

    logging.info('End of Velvet Assembly program')
