#!/usr/bin/env python3

from Bio import SeqIO
import csv
import sys

# set the field size limit to sys maximum
csv.field_size_limit(sys.maxsize)

# snakemake files
names_file = snakemake.params['names']
fasta_file = snakemake.input[0]
outfile = snakemake.output[0]

# generate a dict of read name to new read name
with open(names_file) as csvfile:
    csvreader = csv.reader(csvfile, delimiter='\t')
    header_to_newname = {x[0]: '{0}_{1}'.format(
        x[0],
        len(x[1].split(","))) for x in csvreader}

# read the fasta
handle = SeqIO.parse(fasta_file, 'fasta')
records = [x for x in handle]

# reheader
for rec in records:
    rec.id = header_to_newname[rec.id]
    rec.name = ''
    rec.description = ''

# write output
SeqIO.write(records, outfile, 'fasta')
