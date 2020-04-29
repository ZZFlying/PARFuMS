#!/usr/bin/env python3
from shutil import rmtree
from warnings import warn
from tempfile import TemporaryDirectory


class AutoTempDir(TemporaryDirectory):

    def __init__(self, suffix=None, prefix=None, dir=None, auto_del=True):
        TemporaryDirectory.__init__(self, suffix, prefix, dir)
        self._del = auto_del

    @classmethod
    def _cleanup(cls, name, warn_message):
        rmtree(name)
        warn(warn_message, ResourceWarning)

    def cleanup(self):
        if self._finalizer.detach() and self._del:
            rmtree(self.name)
