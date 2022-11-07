# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

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
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('sequence')
    parser.add_argument('--reverse', action='store_true')
    return parser.parse_args()

def main() -> None:
    args = parse()
    inp = args.sequence
    print('Input')
    print(f"5'{' ' * (len(inp)-3)}3'")
    print(inp + '\n')
    c = complement(inp)
    if args.reverse:
        print(f"5'{' ' * (len(c)-3)}3'")
        print(c[::-1])
    else:
        print(f"3'{' ' * (len(c)-3)}5'")
        print(c)

if __name__ == '__main__':
    main()