#!/usr/bin/env python3

from yaml import safe_load


class Config:
    def __init__(self, config: dict):
        self.fw_file = config['FW_file']
        self.rc_file = config['RC_file']
        self.bc_file = config['BC_file']
        self.work_dir = config['work_dir']
        self.primer_file = config['primer_file']
        self.vector_file = config['vector_file']


def parse_yaml(yaml) -> Config:
    with open(yaml) as steam:
        return Config(safe_load(steam))


if __name__ == '__main__':
    parse_yaml('/home/ubuntu/PARFuMS/examples/configs.yml')
