#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    if len(sys.argv) != 3:
        print('Usage: make_count_matrix.py <quant_dir> <output_tsv>', file=sys.stderr)
        return 2

    quant_dir = Path(sys.argv[1])
    output_tsv = Path(sys.argv[2])
    quant_files = sorted(quant_dir.glob('*.quant.sf'))

    if not quant_files:
        print(f'No quant.sf files found in {quant_dir}', file=sys.stderr)
        return 1

    columns = []
    for quant_file in quant_files:
        sample_id = quant_file.name.replace('.quant.sf', '')
        quant = pd.read_csv(quant_file, sep='\t', usecols=['Name', 'NumReads'])
        columns.append(quant.rename(columns={'NumReads': sample_id}).set_index('Name'))

    matrix = pd.concat(columns, axis=1).fillna(0)
    matrix.index.name = 'transcript_id'
    matrix.to_csv(output_tsv, sep='\t')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
