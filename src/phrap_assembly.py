#!/usr/bin/env python3
import logging
import sys
from shutil import rmtree
from tempfile import mkdtemp
from os import path, mkdir, remove

from sub.parfums_subs import submit_array
from util.singleton_config import Config

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


def get_phrap(seq_file, suffix, param):
    logging.info('Inside get phrap cmd function')
    contig_file = dict()
    script = list()
    # phrap组装后生成三份序列文件，合并输出到一个文件
    # 返回组成生成configs文件的绝对路径
    for (ident, file) in seq_file.items():
        if path.exists(file):
            dirname = path.dirname(file)
            out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
            phrap_file = '{0}.contigs {0}.singlets {0}.problems'.format(file)
            # cmd = 'phrap {} {} && cat $file\.contigs $file\.singlets $file\.problems > $OUTFILE'.format(param, file)
            cmd = str.join(' ', ['phrap', param, file, '&&', 'cat', phrap_file, '>', out_file])
            cmd = cmd + '\n'
            logging.debug('cmd => ' + cmd)
            contig_file[ident] = out_file
            script.append(cmd)
    return script, contig_file


def get_cd_hit_est(contig_file, suffix):
    script = list()
    cdhit_file = dict()
    # 对传入序列文件去重，返回去重文件的绝对路径
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
    # 生成seq文件映射到去重后序列的映射关系文件
    for (ident, mapping_file) in seq_file.items():
        dirname = path.dirname(mapping_file)
        out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
        cmd = 'fr-hit -a {} -d {} -m 30 -o {}\n'.format(mapping_file, cdhit_file[ident], out_file)
        logging.debug('cmd => ' + cmd)
        script.append(cmd)
        frhit_file[ident] = out_file
    return script, frhit_file


def get_link_contigs(cdhit, frhit, temp_dir, suffix):
    script = list()
    phrap_file = dict()
    perl_script = sys.path[0]
    perl_script = path.join(perl_script, 'sub', 'fr-hit_link_light.py')
    for (ident, file) in cdhit.items():
        if path.exists(file):
            dirname = path.dirname(file)
            out_file = path.join(dirname, '{}.{}'.format(ident, suffix))
            cmd = str.join(' ', ['python3', perl_script, frhit[ident], cdhit[ident], temp_dir, ident, out_file])
            cmd += '\n'
            logging.debug('cmd => ' + cmd)
            phrap_file[ident] = out_file
            script.append(cmd)
    return script, phrap_file


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
        if run_type == 'phrap':
            logging.info('Writing Phrap script file')
            script, out_file = get_phrap(seq, suffix, params)
        elif run_type == 'cd-hit-est':
            logging.info('Writing CD-HIT script file')
            script, out_file = get_cd_hit_est(contig, suffix)
        elif run_type == 'fr-hit':
            logging.info('Writing FR-HIT script file')
            script, out_file = get_fr_hit(seq, cdhit, suffix)
        elif run_type == 'link-contigs':
            logging.info('Writing link-contigs script file')
            script, out_file = get_link_contigs(cdhit, frhit, temp_dir, suffix)
        else:
            logging.error('Invalid runType:{}'.format(run_type))

        if len(script) == 0:
            logging.info('Nothing in the command list to submit on cluster!')
        else:
            out.writelines(script)
    return job_file, out_file


def clean_files(files):
    for ident, file in files.items():
        remove(file)


def round_1(work_dir, seq_file, phrap_file):
    logging.info('Phrap Assembly Round-1 Started')
    tempdir = mkdtemp(dir=path.join(work_dir, 'temp'), prefix='PhrapRun1_')

    phrap_param = '-minmatch 25 -maxmatch 40 -bandwidth 1 -minscore 30'
    script, contig_file = make_script('phrap', tempdir, seq=phrap_file, suffix='phrapContigs.fasta', params=phrap_param)
    submit_array(script, 'PhrapRun1', tempdir)

    script, cdhit_file = make_script('cd-hit-est', tempdir, contig=contig_file, suffix='phrap.cd-hit1.fasta')
    submit_array(script, 'CD-hit1', tempdir)
    logging.info('CD-hit-est Run complete')

    script, frhit_file = make_script('fr-hit', tempdir, seq=seq_file, cdhit=cdhit_file, suffix='phrap.Map1.txt')
    submit_array(script, 'FR-hit1', tempdir)
    logging.info('FR-hit Run complete')

    script, phrap_file = make_script('link-contigs', tempdir, frhit=frhit_file, cdhit=cdhit_file, suffix='ForPhrap2.fasta')
    submit_array(script, 'Link-contigs', tempdir)
    logging.info('Link-contigs Run complete')

    clean_files(frhit_file)
    if Config()['auto_del']:
        rmtree(tempdir)
    logging.info('Round-1 completed')
    return phrap_file


def round_2(work_dir, phrap_file):
    logging.info('Phrap Assembly Round-2 Started')
    tempdir = mkdtemp(dir=path.join(work_dir, 'temp'), prefix='PhrapRun2_')

    phrap_param = '-minmatch 25 -maxmatch 40 -bandwidth 1 -minscore 30'
    script, contig_file = make_script('phrap', tempdir, seq=phrap_file, suffix='phrapContigs2.fasta', params=phrap_param)
    submit_array(script, 'PhrapRun2', tempdir)

    script, cdhit_file = make_script('cd-hit-est', tempdir, contig=contig_file, suffix='phrap.fasta')
    submit_array(script, 'CD-hit2', tempdir)
    logging.info('CD-hit-est Run complete')

    if Config()['auto_del']:
        rmtree(tempdir)
    logging.info('Round-2 completed')


def main(work_dir, idents):
    logging.info('Step 5: PHRAP ASSEMBLY STARTED')
    seq_file = read_inputSequences(work_dir, idents, 'noVector.fasta')
    phrap_file = read_inputSequences(work_dir, idents, 'ForPhrap1.fasta')
    phrap_file = round_1(work_dir, seq_file, phrap_file)
    round_2(work_dir, phrap_file)
