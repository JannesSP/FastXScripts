# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from Bio import SeqIO
import os

ACCURATE = 'ACGTU'
AMBIGUOUS = 'KMRSWYN'

IUPAC_DNA = {
    'A':'A',
    'C':'C',
    'G':'G',
    'T':'T',
    'N':'ACGT',
    'Y':'CT',
    'R':'AG',
    'S':'GC',
    'W':'AT',
    'M':'AC',
    'K':'GT'}

IUPAC_RNA = {
    key.replace('T', 'U') : item.replace('T', 'U') for key, item in IUPAC_DNA.items()
}

IUPAC = None

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description='What the fasta will analyse your reference fasta sequence'
    )
    parser.add_argument('FASTA_or_SEQ', type=str, help='FASTA reference file or sequence')
    parser.add_argument('--rna', action='store_true', help='switch to RNA if reference FASTA contains RNA')
    return parser.parse_args()

def count_bases(fasta_sequence : str) -> dict:
    '''
    Counts the IUPAC characters in a FASTA sequence and returns the counts.
    Raises a KeyError if base unknown.

    Parameters
    ----------
    fasta_sequence : str
        The reference FASTA sequence as a string

    Returns
    -------
    counts : dict
        Counts of all nucleotide IUPAC characters
    '''
    counts = {symbol:0 for symbol in IUPAC_DNA}
    for base in fasta_sequence:
        counts[base] += 1
    return counts

def get_seq_content(counts : dict) -> dict:
    '''
    Analyses counts for sequence content

    Parameters
    ----------
    counts : dict
        Counts of all nucleotide IUPAC characters
    
    Returns
    -------
    content : dict
        total : Number of counted bases
        accurate : Number of counted A, C, G or T/U
        ambiguous : Number of counted ambiguous bases like 'N'
        AT : Number counted A and T bases
        GC : Number of counted G and C bases
    '''
    acc = 0 # accurate counts
    amb = 0 # ambiguous counts
    at = 0
    gc = 0

    for character in counts:
    
        if character in ACCURATE:
            acc += counts[character]
        else:
            amb += counts[character]

        if character in 'ATW':
            at += 1
        elif character in 'GCS':
            gc += 1

    return {'total':amb+acc, 'accurate':acc, 'ambiguous':amb, 'AT':at, 'GC':gc}

def output(content : dict, id : str = None) -> None:
    if id is not None:
        print(id)
    print(content)

def main() -> None:
    args = parse()
    fasta = args.fasta
    rna = args.rna

    global IUPAC
    if rna:
        IUPAC = IUPAC_RNA
    else:
        IUPAC = IUPAC_DNA

    # provided FASTA file
    if fasta.endswith('.fa') or fasta.endswith('.fasta'):
        assert os.path.exists(fasta) and os.path.isfile(fasta)

        for record in SeqIO.parse(fasta, format):
            output(get_seq_content(count_bases(str(record.seq))), record.id)

    # provided sequence
    else:
        output(get_seq_content(count_bases(fasta)))

if __name__ == '__main__':
    main()