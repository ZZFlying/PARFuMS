#!/usr/bin/env python3


def load_barcode(file):
    try:
        with open(file) as IN:
            # 尝试读取样本的barcode和标识
            # 跳过#号开头的行数
            barcodes = dict()
            for line in IN.readlines():
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                barcode, ident = line.split()
                barcode = barcode.upper()
                length = len(barcode)
                barcodes[barcode] = ident
            idents = list(barcodes.values())
            return barcodes, idents, length
    except FileNotFoundError:
        pass
