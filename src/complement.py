# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import os
from Bio import SeqIO
from Bio.Seq import Seq

COMPLEMENT = {
    'A':'T',
    'C':'G',
    'G':'C',
    'T':'A',
    'N':'N',
    'Y':'R',
    'R':'Y',
    'S':'S',
    'W':'W',
    'M':'K', 
    'K':'M'
    }

def complement(seq):
    ret = ''
    for b in seq:
        try:
            ret += COMPLEMENT.get(b)
        except KeyError:
            print(f'Base {b} unknown, no complement found!')
            exit(1)
    return ret

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help='Translate nucleotide sequences to its complement'
    )
    parser.add_argument('sequences', help='Input sequence separated with "," or fasta file')
    parser.add_argument('--reverse', action='store_true', help='Use to print 3\'->5\' sequence.')
    return parser.parse_args()

def main() -> None:
    args = parse()
    inp : str = args.sequences
    rev : bool = args.reverse

    if inp.endswith('.fa') or inp.endswith('.fasta'):
        assert os.path.exists(inp) and os.path.isfile(inp)
        outfile = os.path.join(f'{os.path.splitext(inp)[0]}_reverse-complement{os.path.splitext(inp)[1]}') if rev else os.path.join(f'{os.path.splitext(inp)[0]}_complement{os.path.splitext(inp)[1]}')
        records = []

        for rec in SeqIO.parse(inp, 'fasta'):
            c = complement(rec.seq)
            rec.seq = Seq(c[::-1]) if rev else Seq(c)
            rec.id += '_reverse-complement' if rev else '_complement'
            records.append(rec)
            rec.description = ''

        SeqIO.write(records, outfile, 'fasta')

    else:
        for seq in inp.split(','):

            print('Input')
            print(f"5'{' ' * (len(seq)-3)}3'")
            print(seq + '\n')

            c = complement(seq)
            if rev:
                print(f"5'{' ' * (len(c)-3)}3'")
                c = c[::-1]
            else:
                print(f"3'{' ' * (len(c)-3)}5'")
            print(c)

            if c == seq:
                print('PALINDROM!')

if __name__ == '__main__':
    main()