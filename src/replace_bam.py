#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

'''
After using `replace_fastx.py` a list of replaced base positions per sequence is stored.
`replace_bam.py` now fixes the altered mapped reads.
'''

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import pysam

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("inbam")
    parser.add_argument("outbam")
    parser.add_argument("readreplaceCSV")
    return parser.parse_args()

def main() -> None:
    args = parse()
    inbam = pysam.AlignmentFile(args.inbam, 'rb')
    outbam = pysam.AlignmentFile(args.outbam, 'wb', template=inbam)
    readreplaceDF = pd.read_csv(args.readreplaceCSV)

    for ridx, read in enumerate(inbam):
        if (ridx + 1) % 10 == 0:
            print(f'Read {ridx+1}', end='\r')
        readreplacedBases = readreplaceDF[readreplaceDF['readid'] == read.query_name]        
        qualities = read.query_qualities
        query_sequence = list(read.query_sequence)

        for _, replacedPos in readreplacedBases.iterrows():
            query_sequence[replacedPos['position']] = replacedPos['sourcebase']

        read.query_sequence = ''.join(query_sequence)
        read.query_qualities = qualities
        outbam.write(read)

    print(f'Read {ridx+1}')

if __name__ == '__main__':
    main()