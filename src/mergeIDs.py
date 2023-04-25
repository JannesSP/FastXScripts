#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace, FileType
from os.path import exists
from io import TextIOWrapper

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='Write IDs to output file that is present in every provided ID file.'
    )
    parser.add_argument('file', type=FileType('r'), nargs='+', help='Provide a list of ID files')
    parser.add_argument('method', choices=['intersect', 'union'], default='intersect', help='Merge method for lists of IDs')
    parser.add_argument('-o', '--outfile', type=str, help='File to write output IDs', required=True)
    parser.add_argument('-i', '--ignore', action='store_true', default=False)
    return parser.parse_args()

def main() -> None:
    args = parse()
    # This is already a list of TextIOWrappers
    files = args.file
    method = args.method
    outfile = args.outfile
    if not args.ignore:
        assert not exists(outfile), f'{outfile} already exists!'
    assert len(files) > 0
    if method == 'intersect':
        ids = intersect(files)
    elif method == 'union':
        ids = union(files)
    else:
        print('Unknown Merging Method!')
        exit(1)
    write(ids, outfile)

def write(ids : set, outfile : str) -> None:
    '''
    Write set of IDs to given file.

    Parameters
    ----------
    ids: set
    outfile : str
    '''
    with open(outfile, 'w') as w:
        for id in ids:
            w.write(f'{id}\n')

def intersect(files : list) -> set:
    '''
    Iterates through files and returns the set of IDs that are present in all files simultaneously.

    Parameters
    ----------
    files : list

    Returns
    -------
    IDs : set
    '''
    for i, file in enumerate(files):
        tset = set()
        # files is already a list of TextIOWrappers, safety check for test.py
        if file is not TextIOWrapper:
            file = open(file, 'r')
        for line in file:
            line : str = line.strip()
            tset.add(line)
        if file is TextIOWrapper:
            file.close()
        iset = tset.copy() if i == 0 else iset.intersection(tset)

    return iset

def union(files : list) -> set:
    '''
    Iterates through files and returns the union set of IDs from all files.

    Parameters
    ----------
    files : list

    Returns
    -------
    IDs : set
    '''
    uset = set()
    for file in files:
        # files is already a list of TextIOWrappers, safety check for test.py
        if file is not TextIOWrapper:
            file = open(file, 'r')
        for line in file:
            line : str = line.strip()
            uset.add(line)
        if file is TextIOWrapper:
            file.close()
    return uset

if __name__ == '__main__':
    main()