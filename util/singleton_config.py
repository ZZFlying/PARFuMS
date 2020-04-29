#!/usr/bin/env python3

from yaml import safe_load


def Singleton(cls):
    _instance = {}

    def _singleton(*args, **kwargs):
        if cls not in _instance:
            _instance[cls] = cls(*args, **kwargs)
        return _instance[cls]

    return _singleton


def read_config(config_file):
    with open(config_file) as config_in:
        config = safe_load(config_in)
        return config


# 装饰器实现单例模式
@Singleton
class Config(object):
    def __init__(self, config_file=None):
        self._config = read_config(config_file)

    def __getitem__(self, item):
        return self._config[item]

    def __setitem__(self, key, value):
        self._config[key] = value

    def update(self, configs):
        self._config.update(configs)
