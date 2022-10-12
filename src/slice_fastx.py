# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from io import TextIOWrapper
import os
from Bio import SeqIO

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('inFastx', type=str, help='Fastx file from which to slice subsequences')
    parser.add_argument('outFastx', type=str, help='Fastx file to write slices')
    parser.add_argument('--append', action='store_true', help='Appends slices to existing outFastx')
    parser.add_argument('--lowerbound', default=None, type=int, help='Lower bound for slicing area (1-based)')
    parser.add_argument('--upperbound', default=None, type=int, help='Upper bound for slicing area (1-based)')
    parser.add_argument('--position', default=None, type=int, help='Position which to slice (1-based)')
    parser.add_argument('--range', default=None, type=int, help='Range which to slice up- and downstream from the position')
    parser.add_argument('--id', default=None, type=str, help='Fastx ID filter to slice from specific sequence (only works for one ID)')
    return parser.parse_args()

def getSlice(position : int, r : int, lowerbound : int, upperbound : int) -> tuple:
    '''
    Return slice interval for Fastx sequence

    Parameters
    ----------
    position : int
        center position of slice interval
    r : int
        range to include into slice interval
    lowerbound : int
        lower bound of interval that is included
    upperbound : int
        upper bound of interval that is included

    Returns
    -------
    slice : tuple
        interval 0-based [included, excluded)
    '''

    assert (position is None) != (lowerbound is None or upperbound is None)

    # slice python style 0-based [included, excluded)
    slice = [None, None]
    if position is not None:
        if r is not None:
            slice[0] = position - 1 - r
            slice[1] = position + r
        else:
            slice[0] = position - 1
            slice[1] = position
    else:
        if lowerbound is not None:
            slice[0] = lowerbound
        if upperbound is not None:
            slice[1] = upperbound + 1

    assert slice[0] >= 0, f'Lowerbound of slice {slice} is below zero! {slice[0]}'

    return tuple(slice)

def sliceFastx(inFastx : TextIOWrapper, outFastx : TextIOWrapper, slice : tuple, format : str, id : str = None) -> list:
    '''
    Slice sequences and write new Fastx

    Parameters
    ----------
    inFastx : TextIOWrapper
        Fastx handle for incoming sequences
    outFastx : TextIOWrapper
        Fastx handle for outgoing sliced sequences
    slice : tuple
        Tuple containing the 0-based [included, excluded) slice interval
    id : str = None
        Only slice one specific sequence from incoming Fastx file
    format : str = 'fasta'
        format of read and written file, (fasta or fastq)
    '''
    slicedRecords = []

    for record in SeqIO.parse(inFastx, format):
        # with id filter
        if id and record.id == id:
            sliceRecord(record, slice, format)
            slicedRecords.append(record)

        elif not id:
            sliceRecord(record, slice)
            slicedRecords.append(record, format)

    SeqIO.write(slicedRecords, outFastx, format)
    return slicedRecords

def sliceRecord(record : SeqIO.SeqRecord, slice : tuple, format : str) -> None:
    assert len(record.seq) >= slice[1], f'Slice {slice} too large for sequence {record.id} with length {len(record.seq)}'
    record.description += f' sliced=({slice[0]+1},{slice[1]})'
    if format == 'fastq':
        phred_quality = record.letter_annotations['phred_quality'][slice[0] : slice[1]]
        del record.letter_annotations['phred_quality']
    record.seq = record.seq[slice[0] : slice[1]]
    if format == 'fastq':
        record.letter_annotations['phred_quality'] = phred_quality


def main() -> None:
    args = parse()

    slice = getSlice(args.position, args.range, args.lowerbound, args.upperbound)
    id = args.id
    append = args.append

    assert append or os.path.exists(args.outFastx), f'{args.outFastx} already exists! Use a different name or --append'

    with open(args.inFastx, 'r') as inFastx:
        with open(args.outFastx, 'a' if append else 'w') as outFastx:
            sliceFastx(inFastx, outFastx, slice, id)

if __name__ == '__main__':
    main()