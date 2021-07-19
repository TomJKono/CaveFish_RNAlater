#!/usr/bin/env python
"""MicroFinder from SEM, with some modifications to report a 1/0 for each gene
in a supplied FASTA, for whether the gene has a >=6bp homopolymer repeat, or
whether the gene has a longer more complex repeat. Takes one argument:
    1) FASTA file.
Requires Biopython."""

import sys
import re
try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write('This script requries Biopython.\n')
    exit(1)


def classify_repeat(r):
    """Classify a repeat object based on its span and the complexity of the
    match. Will return one of 'H' (homopolymer repeat), 'S' (simple sequence
    repeat), or 'N' (none)."""
    if not r:
        return 'N'
    coord = r.span()
    sp = coord[1] - coord[0]
    unit = r.groups()[0]
    if len(unit) == 1 and sp >= 6:
        return 'H'
    elif len(unit) > 1 and sp >= 3*len(unit):
        return 'S'
    else:
        return 'N'


def main(fasta):
    """Main function."""
    # This is the regular expression that we are searching for
    re_search = re.compile(r'([ATCG]+?)\1{1,}')
    # Print a header
    print('gene_id SSR HPR')
    # Open a path to the FASTA file
    seqhandle = open(fasta, 'r')
    for fastaseq in SeqIO.parse(seqhandle, 'fasta'):
        s = str(fastaseq.seq)
        gid = fastaseq.id.split('|')[0]
        repeats = re.finditer(re_search, s)
        has_ssr = '0'
        has_hom = '0'
        for r in repeats:
            r_type = classify_repeat(r)
            if r_type == 'H':
                has_hom = '1'
            if r_type == 'S':
                has_ssr = '1'
            if has_hom == '1' and has_ssr == '1':
                break
        print(gid, has_ssr, has_hom)
        sys.stderr.write(gid + '\n')
    return


if len(sys.argv) != 2:
    print("""MicroFinder from SEM, with some modifications to report a 1/0 for each gene
in a supplied FASTA, for whether the gene has a >=6bp homopolymer repeat, or
whether the gene has a longer more complex repeat. Takes one argument:
    1) FASTA file.
Requires Biopython.""")
    exit(1)
else:
    main(sys.argv[1])
