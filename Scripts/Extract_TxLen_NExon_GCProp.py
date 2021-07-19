#!/usr/bin/evn python
"""Extract the length, number of exons, and the GC content of the transcripts
for all genes in the cave fish v1.02 genome annotation. This script parses
a FASTA from Ensembl Biomart with the following header information:
    - Sequences: cDNA sequences
    - Gene stable ID
    - Transcript stable ID
    - Exon stable ID
Requires Biopython. Takes one argument:
    1) Biomart export FASTA (gzipped)"""

import sys
import gzip
import pprint
try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write('This script requires Biopython.\n')
    exit(1)

# Make a dictionary to hold the gene-based data
gene_dict = {}
with gzip.open(sys.argv[1], 'rt') as f:
    for rec in SeqIO.parse(f, 'fasta'):
        gene_id, tx_id, exon_list = rec.id.split('|')
        # Count how many exons
        n_exon = len(exon_list.split(';'))
        tx_len = len(rec.seq)
        # GC content = (G+C)/(A+T+C+G)
        gc_prop = (rec.seq.count('G') + rec.seq.count('C')) / (rec.seq.count('G') + rec.seq.count('C') + rec.seq.count('A') + rec.seq.count('T'))
        if gene_id not in gene_dict:
            gene_dict[gene_id] = [(tx_id, tx_len, n_exon, gc_prop)]
        else:
            gene_dict[gene_id].append((tx_id, tx_len, n_exon, gc_prop))

# Then, go through the gene dict and sort it on transcript length, to take the
# longest transcript for each gene
longest_gene_dict = {}
for gene in gene_dict:
    k = sorted(gene_dict[gene], key=lambda x: x[1], reverse=True)
    longest_gene_dict[gene] = k[0]

# Print out a header
print(','.join(['GeneID', 'Longest.TxID', 'Longest.TxLen', 'Longest.NExon', 'Longest.GCProp']))
for gene in sorted(longest_gene_dict):
    toprint = [gene] + [str(i) for i in longest_gene_dict[gene]]
    print(','.join(toprint))
