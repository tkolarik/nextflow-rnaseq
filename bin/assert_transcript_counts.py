#!/usr/bin/env python3

import csv
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) < 3:
        print('Usage: assert_transcript_counts.py <counts_tsv> <expected_sample> [<expected_sample> ...]', file=sys.stderr)
        return 2

    counts_tsv = Path(sys.argv[1])
    expected_samples = sys.argv[2:]

    if not counts_tsv.exists():
        print(f'Missing count matrix: {counts_tsv}', file=sys.stderr)
        return 1

    with counts_tsv.open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        try:
            header = next(reader)
        except StopIteration:
            print(f'Empty count matrix: {counts_tsv}', file=sys.stderr)
            return 1

        if header[0] != 'transcript_id':
            print(f'Expected first column transcript_id, found {header[0]!r}', file=sys.stderr)
            return 1

        observed_samples = header[1:]
        if observed_samples != expected_samples:
            print(f'Expected samples {expected_samples}, found {observed_samples}', file=sys.stderr)
            return 1

        rows = list(reader)

    if not rows:
        print('Count matrix has no transcript rows', file=sys.stderr)
        return 1

    nonzero = False
    for row in rows:
        if len(row) != len(header):
            print(f'Row has {len(row)} columns, expected {len(header)}: {row}', file=sys.stderr)
            return 1
        for value in row[1:]:
            count = float(value)
            if count > 0:
                nonzero = True

    if not nonzero:
        print('Count matrix contains no positive counts', file=sys.stderr)
        return 1

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
