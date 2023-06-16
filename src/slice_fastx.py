    # author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from io import TextIOWrapper
import os
from Bio import SeqIO

FORMATS = {
    '.fa' : 'fasta',
    '.fasta' : 'fasta',
    '.fq' : 'fastq',
    '.fastq' : 'fastq',
    }

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('inFastx', type=str, help='Fastx file from which to slice subsequences')
    parser.add_argument('outFastx', type=str, help='Fastx file to write slices')
    parser.add_argument('--append', action='store_true', help='Appends slices to existing outFastx')
    parser.add_argument('--lowerbound', default=None, type=int, help='Lower bound for slicing area (1-based)')
    parser.add_argument('--upperbound', default=None, type=int, help='Upper bound for slicing area (1-based)')
    parser.add_argument('--center', default=None, type=int, help='Center position which to slice (1-based)')
    parser.add_argument('--range', default=None, type=int, help='Range which to slice up- and downstream from the position')
    parser.add_argument('--slice_start', default=None, type=int, help='Slice number of nucleotides from start of reads')
    parser.add_argument('--slice_end', default=None, type=int, help='Slice number of nucleotides from end of reads')
    # TODO read input ids file
    parser.add_argument('--id', default=None, type=str, help='Fastx ID filter to slice from specific sequence (only works for one ID)')
    return parser.parse_args()

def getSliceRegion(position : int, range : int, lowerbound : int, upperbound : int) -> tuple:
    '''
    Return desired interval of Fastx sequence

    Parameters
    ----------
    position : int
        center position of slice interval
    range : int
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
        if range is not None:
            slice[0] = position - 1 - range
            slice[1] = position + range
        else:
            slice[0] = position - 1
            slice[1] = position
    else:
        if lowerbound is not None:
            slice[0] = lowerbound - 1
        if upperbound is not None:
            slice[1] = upperbound

    assert slice[0] >= 0, f'Lowerbound of slice {slice} is below zero! {slice[0]}'

    return tuple(slice)

def sliceFastx(inFastx : str, outFastx : str, slice : tuple, id : str = None) -> list:
    '''
    Slice sequences and write new Fastx

    Parameters
    ----------
    inFastx : str
        Fastx file path for incoming sequences
    outFastx : str
        Fastx file path for outgoing sliced sequences
    slice : tuple
        Tuple containing the 0-based [included, excluded) slice interval
    id : str = None
        Only slice one specific sequence from incoming Fastx file
    '''
    slicedRecords = []
    informat = FORMATS[os.path.splitext(inFastx)[1]].lower()

    for record in SeqIO.parse(inFastx, informat):
        # with id filter
        if id and record.id == id:
            sliceRecord(record, slice, informat)
            slicedRecords.append(record)

        elif not id:
            sliceRecord(record, slice, informat)
            slicedRecords.append(record)

    SeqIO.write(slicedRecords, outFastx, FORMATS[os.path.splitext(outFastx)[1].lower()])
    return slicedRecords

def sliceRecord(record : SeqIO.SeqRecord, slice : tuple, format : str) -> None:
    assert len(record.seq) >= slice[1], f'Slice {slice} too large for sequence {record.id} with length {len(record.seq)}'
    a = slice[0]
    b = len(record.seq) if slice[1] == -1 else slice[1]
    record.description += f' sliced=({len(record.seq) + a + 1 if a < 0 else a+1},{b})'
    if format == 'fastq':
        phred_quality = record.letter_annotations['phred_quality'][a : b]
        del record.letter_annotations['phred_quality']
    record.seq = record.seq[a : b]
    if format == 'fastq':
        record.letter_annotations['phred_quality'] = phred_quality

def slice_start(num_of_bases : int) -> tuple:
    '''
    Slice sequences and write new Fastx
    
    Parameters
    ----------
    num_of_bases : int
    
    Returns
    -------
    slice : tuple
        interval (0, num_of_bases) 0-based for slice interval [0, num_of_bases)
    '''
    return (0, num_of_bases)

def slice_end(num_of_bases : int) -> tuple:
    '''
    Slice sequences and write new Fastx
    
    Parameters
    ----------
    num_of_bases : int
    
    Returns
    -------
    slice : tuple
        interval (num_of_bases, -1) 0-based for slice interval [len(read) - num_of_bases, len(read))
    '''
    return (-num_of_bases, -1)

def main() -> None:
    args = parse()
    id = args.id
    append = args.append
    assert append or not os.path.exists(args.outFastx), f'{args.outFastx} already exists! Use a different name or --append'
    assert os.path.splitext(args.inFastx)[1].lower() in FORMATS, 'Unknown format of input file'
    assert os.path.splitext(args.outFastx)[1].lower() in FORMATS, 'Unknown format of output file'

    if args.slice_start is not None:
        slice = slice_start(args.slice_start)
    elif args.slice_end is not None:
        slice = slice_end(args.slice_end)
    else:
        slice = getSliceRegion(args.position, args.range, args.lowerbound, args.upperbound)

    sliceFastx(args.inFastx, args.outFastx, slice, id)

if __name__ == '__main__':
    main()
