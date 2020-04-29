#!/usr/bin/env python3
import logging

from os import system
from concurrent.futures import ThreadPoolExecutor, as_completed

from src.singleton_config import Config


def compress(filename):
    cmd = 'gzip -f {}'.format(filename)
    system(cmd)
    return filename


def compress_multi(files_list):
    threads = Config()['thread']
    is_gzip = Config()['gzip']
    if is_gzip:
        pool = ThreadPoolExecutor(threads)
        result = list()
        for file in files_list:
            logging.info('Gzip  {}'.format(file))
            result.append(pool.submit(compress, file))
        for future in as_completed(result):
            gz_file = future.result()
            gz_file = gz_file + '.gz'
            logging.info(gz_file + ' had created.')
