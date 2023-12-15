# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

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
        description='What the alignment will analyse your pairwise sequence alignment'
    )
    parser.add_argument('aln', type=str, help='Pairwise sequence alignment')
    parser.add_argument('--rna', action='store_true', help='switch to RNA if reference FASTA contains RNA')
    return parser.parse_args()

def output(content : dict) -> None:
    for key, item in content.items():
        print(f'{key}: {item}, ', end='')
    muts = content["muts"] + content["dels"] + content["ins"]
    print(f'Alignmentsize: {content["size"]}, Identity: {1-(muts/content["size"])}, Substitutions: {content["muts"]}, Insetions: {content["ins"]}, Deletions: {content["dels"]}')

def compare_pair(recs : dict) -> dict:
    rec_1, rec_2 = recs.keys()
    muts=dels=ins=0
    size=len(recs[rec_1])
    for b1, b2 in zip(recs[rec_1], recs[rec_2]):
        if b1 == '-':
            dels+=1
        elif b2 == '-':
            ins+=1
        elif b1 != b2:
            muts+=1
    return {'muts':muts, 'dels':dels, 'ins':ins, 'size':size}

def compare_multi(recs : dict) -> dict:
    size = len(recs[list(recs.keys())[0]])
    df = pd.DataFrame(columns=['position', 'type'])
    muts = indels = 0
    
    for pos in range(size):
        type = ''
        bases = []

        for id in recs:
            seq = recs[id]
            bases.append(seq[pos])

        bases = np.unique(bases)

        if len(bases) == 1:
            continue

        if '-' in bases:
            type='gap'
        
        bases = np.delete(bases, np.where(bases == '-'))
        
        if len(bases):
            if type == 'gap': 
                type='both'
            else:
                type='substitution'

        new_entry = pd.DataFrame({
            'position' : [pos],
            'type' : [type],
        })
        df = pd.concat((new_entry, df), ignore_index=True)
    
    plot_multi(df)
    return {'muts':muts, 'dels':indels, 'ins':indels, 'size':size}

def plot_multi(df : pd.DataFrame) -> None:
    sns.stripplot(data = df, x = 'position', hue='type')
    plt.tight_layout()
    plt.savefig('alignment_muts.pdf')
    plt.savefig('alignment_muts.png')

def main() -> None:
    args = parse()
    aln = args.aln
    rna = args.rna

    global IUPAC
    if rna:
        IUPAC = IUPAC_RNA
    else:
        IUPAC = IUPAC_DNA

    recs = {}
    for record in SeqIO.parse(aln, 'fasta'):
        recs[record.id] = str(record.seq)

    if len(recs) == 2:
        output(compare_pair(recs))
    elif len(recs) > 2:
        output(compare_multi(recs))
    else:
        print("ERROR not enough sequences found")
        exit(1)

if __name__ == '__main__':
    main()
