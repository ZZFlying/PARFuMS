#!/usr/bin/env python3


def load_barcode(file):
    try:
        with open(file) as IN:
            barcodes = dict()
            for line in IN.readlines():
                if line.startswith('#'):
                    continue
                barcode, ident = line.split()
                barcode = barcode.upper()
                length = len(barcode)
                barcodes[barcode] = ident
            idents = list(barcodes.values())
            return barcodes, idents, length
    except FileNotFoundError:
        pass
