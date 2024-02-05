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

def revComp(seq : str) -> str:
    complement_dict = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    rev_seq = seq[::-1]
    rev_comp_seq = ''.join(complement_dict[base] for base in rev_seq)
    return rev_comp_seq

def main() -> None:
    args = parse()
    inbam = pysam.AlignmentFile(args.inbam, 'rb')
    outbam = pysam.AlignmentFile(args.outbam, 'wb', template=inbam)
    readreplaceDF = pd.read_csv(args.readreplaceCSV)

    for ridx, read in enumerate(inbam):
        if (ridx + 1) % 10 == 0:
            print(f'Read {ridx+1}', end='\r')
        
        if read.is_supplementary:
            # TODO how to handle these?
            continue

        readreplacedBases = readreplaceDF[readreplaceDF['readid'] == read.query_name]
        qualities = read.query_qualities
        query_sequence = list(str(read.get_forward_sequence()))

        try:
            for _, replacedPos in readreplacedBases.iterrows():
                assert query_sequence[replacedPos['position']] == replacedPos['targetbase']
                query_sequence[replacedPos['position']] = replacedPos['sourcebase']
        except IndexError as e:
            print("IndexError:")
            print(read.query_name, len(query_sequence), read.query_alignment_start)
            print(read.get_forward_sequence())
            print(replacedPos['position'], query_sequence[replacedPos['position']], replacedPos['targetbase'])
            print(e.with_traceback())
            exit(1)
        except AssertionError as e:
            print("AssertionError:")
            print(read.query_name, len(query_sequence), read.query_alignment_start)
            print(read.get_forward_sequence())
            print(replacedPos['position'], query_sequence[replacedPos['position']], replacedPos['targetbase'])
            print(e.with_traceback())
            exit(1)

        if read.is_reverse:
            read.query_sequence = revComp(''.join(query_sequence))
        else:
            read.query_sequence = ''.join(query_sequence)
        read.query_qualities = qualities
        outbam.write(read)

    print(f'Read {ridx+1}')

if __name__ == '__main__':
    main()