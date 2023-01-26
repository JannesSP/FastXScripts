# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from io import TextIOWrapper
from Bio import SeqIO
import numpy as np

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Filter FASTA or FASTQ file for ids or length of reads'
    )

    parser.add_argument('inFASTX', metavar='FASTX', type=str, help='Multi FASTQ or FASTA file')
    mode = parser.add_mutually_exclusive_group(required = True)
    mode.add_argument('-i', '--read_ids', metavar='IDS', type=str, default=None, help='One read ID per line in file, line separated read IDs')
    mode.add_argument('-l', '--long', metavar='LENGTH', type=int, default=None, help='Filter FASTA or FASTQ file for reads given length or longer')
    mode.add_argument('-s', '--short', metavar='LENGTH', type=int, default=None, help='Filter FASTA or FASTQ file for reads given length or shorter')
    parser.add_argument('-o', '--outFASTX', metavar='FASTX', type=str, default = None, help='FASTQ or FASTA file containing provided reads', required=True)
    nt_mode = parser.add_mutually_exclusive_group()
    nt_mode.add_argument('--dna', action='store_true', default=False, help='Convert output sequences to dna (ACGT)')
    nt_mode.add_argument('--rna', action='store_true', default=False, help='Convert output sequences to rna (ACGU)')

    return parser.parse_args()

def main() -> None:
    args = parse()

    inFX = args.inFASTX
    ids = args.read_ids
    outFX = args.outFASTX
    long = args.long
    short = args.short
    dna=args.dna
    rna=args.rna

    print('Filtering', inFX)

    if outFX.lower().endswith('.fa') or outFX.lower().endswith('.fasta'):
        format = 'fasta'
    elif outFX.lower().endswith('.fq') or outFX.lower().endswith('.fastq'):
        format = 'fastq'
    else:
        print('Unknown file format for', outFX)
        print('Must be .fa/.fasta or .fq/.fastq')
        exit(1)

    if ids is not None:
        filterIDs(open(inFX, 'r'), open(ids, 'r'), open(outFX, 'w+'))

    elif long is not None:
        records, longest, shortest = filterLength(inFX, long, 'long')

    elif short is not None:
        records, longest, shortest = filterLength(inFX, short, 'short')

    for record in records:
        if dna:
            record.seq = record.seq.replace('U', 'T')
        if rna:
            record.seq = record.seq.replace('T', 'U')

    SeqIO.write(records, outFX, format)
    print('longest', longest, 'shortest', shortest)

def filterLength(inFX : str, threshold : int, mode : str) -> None:
    func = {
        'long':lambda length, threshold: True if length >= threshold else False,
        'short':lambda length, threshold: True if length <= threshold else False
        }[mode]
    
    if inFX.lower().endswith('.fa') or inFX.lower().endswith('.fasta'):
        format = 'fasta'
    elif inFX.lower().endswith('.fq') or inFX.lower().endswith('.fastq'):
        format = 'fastq'
    else:
        print('Unknown file format for', inFX)
        print('Must be .fa/.fasta or .fq/.fastq')
        exit(1)

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

def filterIDs(inFX : TextIOWrapper, ids : TextIOWrapper, outFX : TextIOWrapper = None) -> tuple:
    '''
    Filters the input FASTA/FASTQ for ids in given list and writes filtered FASTA/FASTQ.
    
    If outFX is not None:
        Writes found reads to outFX file.

    Parameters
    ----------
    inFX : TextIOWrapper
        Input FASTA/FASTQ
    ids : list
        ReadIDs to filter for
    outFX : TextIOWrapper
        Filtered output FASTA/FASTQ

    Returns
    -------
    foundIDs : list
        list of found IDs
    filtered : list
        list of removed IDs
    missingIDs : list
        list of missing IDs
    '''

    ids_list = list(map(lambda id : id[:-1] if '\n' in id else id, ids.readlines()))
    writeRead = False
    foundIDs = []
    filteredIDs = []

    print(f'Looking for {len(ids_list)} ids')

    for idx, line in enumerate(inFX, 1):
        if idx%10000==0:
            print(f'Processing line {idx} ...', end='\r')
        if line.startswith('@') or line.startswith('>'):
            writeRead = False
            readid = line.split()[0][1:]
            if readid in ids_list:
                writeRead = True
                if outFX is not None:
                    outFX.write(line)
                ids_list.remove(readid)
                foundIDs.append(readid)
            else:
                filteredIDs.append(readid)
        else:
            if writeRead and outFX is not None:
                outFX.write(line)

    print(f'Processed line {idx}\t\t')
    print(f'Found {len(foundIDs)} ids')
    print(f'Filtered {len(filteredIDs)} ids')
    print(f'Missed {len(ids_list)} ids')

    return foundIDs, filteredIDs, ids_list

if __name__ == '__main__':
    main()