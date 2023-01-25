# FastXScripts

Scripts that I needed to change Fasta/FastQ files during my PhD

## Usage

### filter_fastx.py

Filter Fasta/FastQ file for IDs or length of read

```{r}
usage: filter_fastx.py [-h] (-i IDS | -l LENGTH | -s LENGTH) [-o FASTX] FASTX

Filter FASTA or FASTQ file for ids or length of reads

positional arguments:
  FASTX                 Multi FASTQ or FASTA file

optional arguments:
  -h, --help            show this help message and exit
  -i IDS, --read_ids IDS
                        One read ID per line in file, line separated read IDs (default: None)
  -l LENGTH, --long LENGTH
                        Filter FASTA or FASTQ file for reads given length or longer (default: None)
  -s LENGTH, --short LENGTH
                        Filter FASTA or FASTQ file for reads given length or shorter (default: None)
  -o FASTX, --outFASTX FASTX
                        FASTQ or FASTA file containing provided reads (default: None)
```

### slice_fastx.py

Slice subsequences by their position from reads in Fasta/FastQ.

```
usage: slice_fastx.py [-h] [--append] [--lowerbound LOWERBOUND] [--upperbound UPPERBOUND] [--position POSITION] [--range RANGE] [--id ID] inFastx outFastx

positional arguments:
  inFastx               Fastx file from which to slice subsequences
  outFastx              Fastx file to write slices

options:
  -h, --help            show this help message and exit
  --append              Appends slices to existing outFastx (default: False)
  --lowerbound LOWERBOUND
                        Lower bound for slicing area (1-based) (default: None)
  --upperbound UPPERBOUND
                        Upper bound for slicing area (1-based) (default: None)
  --position POSITION   Position which to slice (1-based) (default: None)
  --range RANGE         Range which to slice up- and downstream from the position (default: None)
  --id ID               Fastx ID filter to slice from specific sequence (only works for one ID) (default: None)
```

### complement.py

Translate nucleotide sequences from terminal or fasta files.

```
usage: complement.py [-h] [--reverse] sequences

positional arguments:
  sequences   Input sequence separated with "," or fasta file

optional arguments:
  -h, --help  show this help message and exit
  --reverse   Use to print 3'->5' sequence. (default: False)
  --rna       Translate RNA sequences (default: False)
```

### wtf.py

What the fasta will analyse given sequences for their content like number of bases, the AT and GC content, the number of ambiguous bases (e.g. N).

```
What the fasta will analyse your reference fasta sequence

positional arguments:
  FASTA_or_SEQ  FASTA reference file or sequence

options:
  -h, --help    show this help message and exit
  --rna         switch to RNA if reference FASTA contains RNA (default: False)
```