#!/usr/bin/env python3

import sys
from os import path
from collections import defaultdict


def cross_match_script(split_file, map_file, tempdir, py_script, suffix):
    results = defaultdict(list)
    script = 'cross_match_script.sh'
    script = path.join(tempdir, script)
    sub_folder = path.join(sys.path[0], 'sub')
    py_script = path.join(sub_folder, py_script)
    with open(script, 'w') as output:
        for (ident, files) in split_file.items():
            for file in files:
                cmd = 'cross_match {} {}'.format(file, map_file)
                cmd += ' -gap1_only -minmatch 6 -minscore 10 -gap_init -3'
                cmd += ' | python3 {} {}'.format(py_script, file)
                cmd += '\n'
                output.write(cmd)
                out_file = '{}.{}'.format(file, suffix)
                results[ident].append(out_file)
    return script, results
