#!/usr/bin/env python3
import gzip


class record:
    def __init__(self, id, seq):
        self.id = id
        self._seq = seq

    def __str__(self):
        return self.parse()

    def parse(self, length=50):
        _seq = [self._seq[i:i + length] for i in range(0, len(self._seq), length)]
        self._seq = str.join('\n', _seq) + '\n'
        return self.id + self._seq


class fastaReader:

    def __init__(self, handle):
        _open = gzip.open if handle.split('.')[-1] == 'gz' else open
        self._handle = _open(handle, 'rt')

    def __iter__(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __next__(self):
        try:
            self._read = self._read_one_record()
            return self._read
        except StopIteration:
            raise StopIteration

    def _read_one_record(self):
        _id = self._handle.readline()
        while not _id.startswith('>'):
            if not _id:
                raise StopIteration
            _id = self._handle.readline()

        _seq = ''
        pos = self._handle.tell()
        _sub_seq = self._handle.readline()
        while not _sub_seq.startswith('>'):
            _seq += _sub_seq.strip()
            pos = self._handle.tell()
            _sub_seq = self._handle.readline()
            if not _sub_seq:
                self._handle.seek(pos)
                return record(_id, _seq)
        self._handle.seek(pos)
        return record(_id, _seq)

    def close(self):
        self._handle.close()
