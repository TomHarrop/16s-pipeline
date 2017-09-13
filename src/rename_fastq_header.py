#!/usr/bin/env python3

import gzip
from Bio import SeqIO

fq = snakemake.input[0]
samplename = snakemake.params['samplename']
sampletype = snakemake.params['sampletype']
outfile = snakemake.output[0]

# open read file
with gzip.open(fq, 'rt') as handle:
    fastq = SeqIO.parse(handle, 'fastq-sanger')
    records = [x for x in fastq]
    # rename records
    i = 0
    for rec in records:
        i += 1
        rec.id = '{0}_{1}|{2}'.format(
            samplename, sampletype, str(i))
        rec.name = ''
        rec.description = ''

# write to fasta
SeqIO.write(records, outfile, 'fasta')
