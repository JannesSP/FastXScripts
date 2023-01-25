# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from io import TextIOWrapper
from Bio import SeqIO

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
    parser.add_argument('-o', '--outFASTX', metavar='FASTX', type=str, default = None, help='FASTQ or FASTA file containing provided reads')

    return parser.parse_args()

def main() -> None:
    args = parse()

    inFX = args.inFASTX
    ids = args.read_ids
    outFX = args.outFASTX
    long = args.long
    short = args.short

    if ids is not None:
        filterIDs(open(inFX, 'r'), open(ids, 'r'), open(outFX, 'w+'))

    elif long is not None:
        filterLength(inFX, outFX, long, 'long')

    elif short is not None:
        filterLength(inFX, outFX, short, 'short')

def filterLength(inFX : str, outFX : str, threshold : int, mode : str) -> None:
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
    for i, seq_record in enumerate(infx):
        if (i+1)%100==0:
            print('Checking read', i+1, '\r')
        if func(len(seq_record), threshold):
            out.append(seq_record)
    SeqIO.write(out, outFX, format)

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