# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from src.filter_fastx import filterFX
from src.slice_fastx import sliceFastx, getSlice
import os

testFastq = os.path.join(os.path.dirname(__file__), 'test.fastq')
testFasta = os.path.join(os.path.dirname(__file__), 'test.fasta')
ids = os.path.join(os.path.dirname(__file__), 'ids.txt')

found = set(['4052e08f-635c-419f-acd0-383c7ba40daa', '5e99c5f5-61a9-48e8-9d3c-b74f82e338cc'])
filtered = set(['cf9d662a-850f-49a9-ab55-86ea4e34aa23', '9ca3164f-ce83-4e41-ae86-44a1646aaec4', 'c672e8d3-0b3b-48b5-8170-1ad20e923257'])
missing = set(['c672e8d3-0b3b-48b5-8170-1ad20e923251', 'c6eee8d3-033c-4755-86a0-1ad20e923257'])

def test_filter_fastq():
    foundIDs, filteredIDs, missingIDs = filterFX(open(testFastq, 'r'), open(ids, 'r'), open(os.path.join(os.path.dirname(__file__), 'outfiles', 'fastq_out.fq'), 'w'))
    assert found == set(foundIDs)
    assert filtered == set(filteredIDs)
    assert missing == set(missingIDs)

def test_filter_fasta():
    foundIDs, filteredIDs, missingIDs = filterFX(open(testFasta, 'r'), open(ids, 'r'), open(os.path.join(os.path.dirname(__file__), 'outfiles', 'fastq_out.fa'), 'w'))
    assert found == set(foundIDs)
    assert filtered == set(filteredIDs)
    assert missing == set(missingIDs)

def test_slice_fasta_id():
    records = sliceFastx(open(testFasta, 'r'), open(os.path.join(os.path.dirname(__file__), 'outfiles', 'fasta_out_sliced.fa'), 'w'), (8,15), 'fasta', '4052e08f-635c-419f-acd0-383c7ba40daa')
    assert records[0].seq == 'AGGUAUC'

def test_slice_fastq_id():
    records = sliceFastx(open(testFastq, 'r'), open(os.path.join(os.path.dirname(__file__), 'outfiles', 'fastq_out_sliced.fq'), 'w'), (8,15), 'fastq', '4052e08f-635c-419f-acd0-383c7ba40daa')
    assert records[0].seq == 'AGGUAUC'

# overwrites output file from test_clise_fasta_id()
def test_slice_fasta():
    records = sliceFastx(open(testFasta, 'r'), open(os.path.join(os.path.dirname(__file__), 'outfiles', 'fasta_out_sliced.fa'), 'w'), (8,15), 'fasta')
    seqs = [str(record.seq) for record in records]
    assert seqs == ['AGGUAUC', 'UGAUUUA', 'GUGCCCC', 'ACGUCAC', 'CCCACCC']

def test_getSlice_position_r():
    assert (4, 15) == getSlice(position = 10, r = 5, lowerbound = None, upperbound = None)

def test_getSlice_lower_upper():
    assert (4, 15) == getSlice(position = None, r = None, lowerbound = 5, upperbound = 15)