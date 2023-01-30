#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from io import TextIOWrapper
# import os
# import h5py
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# import scipy
import pysam

from src.filter_fastx import filterLength

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Filter bam file for read lengths and write them to a fasta file, fastq not available yet'
    )
    parser.add_argument('--bam', type=str, help='Mapping bam file', required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    # mode.add_argument('-i', '--read_ids', metavar='IDS', type=str, default=None, help='One read ID per line in file, line separated read IDs')
    mode.add_argument('-l', '--long', metavar='LENGTH', type=int, default=None, help='Filter BAM file for reads given length or longer')
    mode.add_argument('-s', '--short', metavar='LENGTH', type=int, default=None, help='Filter BAM file for reads given length or shorter')
    parser.add_argument('--outfile', type=str, help='fasta file to write reads to', required=True)
    return parser.parse_args()

def main() -> None:
    args = parse()
    bam = args.bam
    long=args.long
    short=args.short
    outfile=args.outfile

    assert outfile.lower().endswith('.fa') or outfile.lower().endswith('.fasta')

    if long is not None:
        filterLength(bam, outfile, long, 'long')
    elif short is not None:
        filterLength(bam, outfile, short, 'short')

def filterLength(bamfile : str, outfile : str, threshold : int, mode : str) -> None:
    func = {
        'long':lambda length: True if length >= threshold else False,
        'short':lambda length: True if length <= threshold else False
        }[mode]

    out = open(outfile, 'w')

    bam = pysam.Samfile(bamfile, "rb")
    for read in bam.fetch():
        if read.is_mapped:
            if func(read.query_alignment_length):
                writeFastx(out, read.query_name, read.query_sequence)
                print(f'Found {read.query_name} with length {read.query_alignment_length}', end = '\r')
    print('Done')

def writeFastx(file : TextIOWrapper, name : str, read : str) -> None:
    file.write(f'>{name}\n{read}\n')

if __name__ == '__main__':
    main()