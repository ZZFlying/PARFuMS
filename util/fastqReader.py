#!/usr/bin/env python3
import gzip


class record:
    def __init__(self, description, seq, plus, phred, id, barcode):
        self._description = description
        self._seq = seq
        self._plus = plus
        self._phred = phred
        self.id = id
        self.barcode = barcode

    def __str__(self):
        return self._description + self._seq + self._plus + self._phred

    def parse(self, toward):
        return self.id + "#0_{}\n".format(toward) + self._seq + self._plus + self._phred


class fastqReader:

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
        _description = self._handle.readline()
        while not _description.startswith('@'):
            if not _description:
                raise StopIteration
            _description = self._handle.readline()
        _id = _description.split()[0]
        _barcode = _description.split(':')[-1].strip()
        _seq = self._handle.readline()
        _plus = self._handle.readline()
        _phred = self._handle.readline()
        return record(_description, _seq, _plus, _phred, _id, _barcode)

    def close(self):
        self._handle.close()
