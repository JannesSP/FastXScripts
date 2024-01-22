#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

'''
This script takes:
- a fasta/fastq file,
- a source base to replace,
- a target base to replace the source base with.

It produces an output fasta/fastq file with the replaced base and a csv file with information about the replaces bases.
'''

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
from Bio import SeqIO
import re

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("fastx", type = str, help = "Input fasta/fastq file.")
    parser.add_argument("srcbase", type = str, help = "Source base to be replaced.")
    parser.add_argument("tgtbase", type = str, help = "Target base that is inserted.")
    parser.add_argument("outdir", type = str, help = "Output directory")
    return parser.parse_args()

def replaceBase(file : str, format : str, srcbase: str, tgtbase : str, outfastx : str, outcsv : str) -> None:

    sequences = []

    with open(outcsv, 'w') as csv:
        csv.write('readid,position,sourcebase,targetbase\n')

        for record in SeqIO.parse(file, format):
            indices = re.finditer(srcbase, str(record.seq))
            record.seq = record.seq.replace(srcbase, tgtbase)
            sequences.append(record)
            outlines = [f'{record.id},{idx.start()},{srcbase},{tgtbase}' for idx in indices]
            if outlines:
                csv.write('\n'.join(outlines) + '\n')
        
    SeqIO.write(sequences, outfastx, format)

def main() -> None:
    args = parse()
    fastx = args.fastx
    srcbase = args.srcbase
    tgtbase = args.tgtbase
    
    if fastx.endswith('fasta') or fastx.endswith('fa') or fastx.endswith('fn'):
        format = 'fasta'
    elif fastx.endswith('fastq') or fastx.endswith('fq'):
        format = 'fastq'
    else:
        print(f'Error: Unknown file extension {os.path.splitext(fastx)[1]}')
        exit(1)

    outfastx = os.path.splitext(fastx)[0] + f'_replaced{srcbase}{tgtbase}.{format}'
    outcsv = os.path.splitext(fastx)[0] + f'_replaced{srcbase}{tgtbase}.csv'
    replaceBase(fastx, format, srcbase, tgtbase, outfastx, outcsv)

if __name__ == '__main__':
    main()