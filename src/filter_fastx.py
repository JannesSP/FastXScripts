# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from io import TextIOWrapper
from Bio import SeqIO
import numpy as np
from os.path import exists

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Filter FASTA or FASTQ file for ids or length of reads'
    )

    parser.add_argument('inFASTX', metavar='inFASTX', type=str, help='Multi FASTQ or FASTA file')
    parser.add_argument('outFASTX', metavar='outFASTX', type=str, help='FASTQ or FASTA file containing provided reads')
    mode = parser.add_mutually_exclusive_group(required = True)
    mode.add_argument('-i', '--read_ids', metavar='IDS', type=str, default=None, help='One read ID per line in file, line separated read IDs')
    mode.add_argument('-l', '--long', metavar='LENGTH', type=int, default=None, help='Filter FASTA or FASTQ file for reads given length or longer')
    mode.add_argument('-s', '--short', metavar='LENGTH', type=int, default=None, help='Filter FASTA or FASTQ file for reads given length or shorter')
    nt_mode = parser.add_mutually_exclusive_group()
    nt_mode.add_argument('--dna', action='store_true', default=False, help='Convert output sequences to dna (ACGT)')
    nt_mode.add_argument('--rna', action='store_true', default=False, help='Convert output sequences to rna (ACGU)')

    return parser.parse_args()

def main() -> None:
    args = parse()

    inFX=args.inFASTX
    ids=args.read_ids
    outFX=args.outFASTX
    long=args.long
    short=args.short
    dna=args.dna
    rna=args.rna

    print('Filtering', inFX)

    # Check Parameters
    if inFX.lower().endswith('.fa') or inFX.lower().endswith('.fasta'):
        informat = 'fasta'
    elif inFX.lower().endswith('.fq') or inFX.lower().endswith('.fastq'):
        informat = 'fastq'
    else:
        print('Unknown input file format for', inFX)
        print('Must be .fa/.fasta or .fq/.fastq')
        exit(1)

    if outFX.lower().endswith('.fa') or outFX.lower().endswith('.fasta'):
        outformat = 'fasta'
    elif outFX.lower().endswith('.fq') or outFX.lower().endswith('.fastq'):
        outformat = 'fastq'
    else:
        print('Unknown output file format for', outFX)
        print('Must be .fa/.fasta or .fq/.fastq')
        exit(2)

    assert not exists(outFX), 'outFASTX already exists!'
    assert exists(inFX), 'inFASTX does not exist!'

    if ids is not None:
        assert exists(ids), 'IDs file does not exist!'
        records, filtered, missed = filterIDs(inFX, informat, open(ids, 'r'))
        print('Found Reads: ', len(records), ', Filtered:, ', len(filtered), ', Unseen IDs: ', len(missed))

    else:
        if long is not None:
            records, longest, shortest = filterLength(inFX, informat, long, 'long')

        elif short is not None:
            records, longest, shortest = filterLength(inFX, informat, short, 'short')
        print('longest', longest, 'shortest', shortest)

    for record in records:
        if dna:
            record.seq = record.seq.replace('U', 'T')
        if rna:
            record.seq = record.seq.replace('T', 'U')

    SeqIO.write(records, outFX, outformat)

def filterLength(inFX : str, format : str, threshold : int, mode : str) -> tuple:
    '''
    Filters the input FASTX for reads with a given length.
    Filtering for shorter or longer reads is determined by the mode.

    Parameters
    ----------
    inFX : str
        input FASTX containing the reads
    threshold : int
        filter for this threshold
    mode : str
        filtering for 'long' or 'short' reads

    Returns
    -------
    out : list
        list of desired reads
    longest : int
        length of longest read in inFX
    shortest : int
        length of shortest read in inFX
    '''
    func = {
        'long':lambda length, threshold: True if length >= threshold else False,
        'short':lambda length, threshold: True if length <= threshold else False
        }[mode]

    infx = SeqIO.parse(inFX, format)
    out = []
    longest = -np.inf
    shortest = np.inf
    for i, seq_record in enumerate(infx):
        if (i+1)%1000==0:
            print('Checking read', i+1, '\tFound', len(out), end='\r')
        if func(len(seq_record), threshold):
            out.append(seq_record)
        longest = max(longest, len(seq_record))
        shortest = min(shortest, len(seq_record))
    print()
    return out, longest, shortest

def filterIDs(inFX : str, format : str, ids : TextIOWrapper) -> tuple:
    '''
    Filters the input FASTX for ids in given list and writes filtered FASTX.
    
    If outFX is not None:
        Writes found reads to outFX file.

    Parameters
    ----------
    inFX : str
        Input FASTA/FASTQ
    format : str
        define the input file format
    ids : TextIOWrapper
        ReadIDs to filter for

    Returns
    -------
    foundRecords : list
        list of found SeqRecords
    removedIDs : list
        list of removed IDs
    missedIDs : list
        list of missed IDs
    '''

    ids_list = set(map(lambda id : id.strip(), ids.readlines()))
    foundRecords = []
    removedIDs = []

    infx = SeqIO.parse(inFX, format)

    print(f'Looking for {len(ids_list)} ids')

    for idx, seq_record in enumerate(infx):
        if (idx+1)%1000==0:
            print(f'Processing line {idx+1} ...', end='\r')

        if seq_record.name in ids_list:
            foundRecords.append(seq_record)
            ids_list.remove(seq_record.name)
        else:
            removedIDs.append(seq_record.name)

    print(f'Processed line {idx}\t\t')
    print(f'Found {len(foundRecords)} ids')
    print(f'Filtered {len(removedIDs)} ids')
    print(f'Missed {len(ids_list)} ids')

    return foundRecords, removedIDs, ids_list

if __name__ == '__main__':
    main()