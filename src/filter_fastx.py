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
    mode.add_argument('-n', '--number', metavar='NUMBER', type=int, default=None, help='Filter FASTA or FASTQ file for given number of reads')
    nt_mode = parser.add_mutually_exclusive_group()
    nt_mode.add_argument('--dna', action='store_true', default=False, help='Convert output sequences to dna (ACGT)')
    nt_mode.add_argument('--rna', action='store_true', default=False, help='Convert output sequences to rna (ACGU)')
    parser.add_argument('-f', '--force', action='store_true', help='Force output overwrite')

    return parser.parse_args()

def main() -> None:
    args = parse()

    inFX=args.inFASTX
    ids=args.read_ids
    outFX=args.outFASTX
    long=args.long
    short=args.short
    number=args.number
    dna=args.dna
    rna=args.rna
    force=args.force

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

    if not force:
        assert not exists(outFX), f'{outFX} already exists!'
    assert exists(inFX), f'{inFX} does not exist!'

    if ids is not None:
        assert exists(ids), f'{ids} file does not exist!'
        records, filtered, missed = filterIDs(inFX, informat, open(ids, 'r'))
        print('Found Reads: ', len(records), ', Filtered:, ', len(filtered), ', Unseen IDs: ', len(missed))

    elif long is not None:
        records, longest, shortest = filterLength(inFX, informat, long, 'long')

    elif short is not None:
        records, longest, shortest = filterLength(inFX, informat, short, 'short')
        print('longest', longest, 'shortest', shortest)
    
    elif number is not None:
        records = filterNum(inFX, informat, number)

    for record in records:
        if dna:
            record.seq = record.seq.replace('U', 'T')
        if rna:
            record.seq = record.seq.replace('T', 'U')

    SeqIO.write(records, outFX, outformat)

def filterNum(inFX : str, format : str, number : int) -> list:
    '''
    Filters (uniformly) randomly drawn reads from given FASTX file.

    Parameters
    ----------
    inFX : str
        input FASTX containing reads
    number : int
        number of randomly chosen reads

    Returns
    -------
    chosen : list
        list with SeqRecords of randomly chosen reads
    '''

    read_dict = SeqIO.to_dict(SeqIO.parse(inFX, format))
    readids = list(read_dict.keys())

    choice = np.random.default_rng().choice(len(readids), size=number, replace=False)
    chosen = []

    for i in choice:
        chosen.append(read_dict[readids[i]])

    return chosen

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

    print(f'Processed lines {idx}\t\t')

    return foundRecords, removedIDs, ids_list

if __name__ == '__main__':
    main()