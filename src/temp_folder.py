#!/usr/bin/env python3
import string as _string
import shutil as _shutil
import weakref as _weakref
import warnings as _warnings
from random import sample as _sample
from os import path as _path, getcwd as _getcwd, mkdir as _mkdir


def mkdtemp(suffix=None, prefix=None, path=None):
    _count = 0
    _candidate = _string.ascii_letters + _string.digits
    _letter_list = [i for i in _candidate]
    while True:
        if _count > 10000:
            raise FileExistsError
        _name = ''.join(_sample(_letter_list, 7))
        if suffix:
            _name = _name + suffix
        if prefix:
            _name = prefix + _name
        if path:
            _dir = _path.join(path, _name)
        else:
            _dir = _path.join(_getcwd(), _name)
        if _path.exists(_dir):
            _count += 1
            continue
        _mkdir(_dir)
        return _dir


class AutoTempdir(object):
    # 生成临时文件夹，根据设置决定失去失去引用后是否删除
    def __init__(self, suffix=None, prefix=None, dir=None, autodel=True):
        self.name = mkdtemp(suffix, prefix, dir)
        self._autodel = autodel
        self._finalizer = _weakref.finalize(
            self, self._cleanup, self.name,
            warn_message="Implicitly cleaning up {!r}".format(self))

    def _cleanup(self, name, warn_message):
        if self._autodel:
            _shutil.rmtree(name)
            _warnings.warn(warn_message, ResourceWarning)

    def __repr__(self):
        return "<{} {!r}>".format(self.__class__.__name__, self.name)
