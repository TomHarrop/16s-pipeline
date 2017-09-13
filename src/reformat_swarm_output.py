#!/usr/bin/env python3

# as far as i can work out, the output of swarm2gut.rb is two files. First, the
# "headers" file just contains a list of records that should be retrieved from
# the dereplicated "unique.fa" file. I think Ruy's script retrieves clusters
# with 10 or more representatives only, prints the header name into a text file
# and then subsets the fasta file on that. The NAMES file should contain a list
# of all the reads in that cluster, generated by parsing the earlier MOTHUR
# output. Because we have biopython I can skip the step of calling
# FaSomeRecords

from Bio import SeqIO
import csv
import re

mothur_names = snakemake.input['mothur_names']
swarm_results = snakemake.input['swarm_results']
dereplicated_fasta = snakemake.input['dereplicated_fasta']
output_names = snakemake.output['names']
output_fasta = snakemake.output['fasta']

# generate a dict of names from the mothur namefile
with open(mothur_names, 'rt') as mothur:
    split_strings = [str(x).rstrip('\n').split('\t') for x in mothur]
    mothur_name_dict = {x[0]: x[1] for x in split_strings}

# read the swarm results
with open(swarm_results, 'rt') as swarm:
    readlists = [str(x).rstrip('\n').split(' ') for x in swarm]

# keep the swarm results with more than 10 representatives
kept_swarms = {x[0]: x for x in readlists
               if len(x) > 9}

# remove the trailing number of reads from the header and the read list
renamed_swarms = {}
for key in kept_swarms:
    new_key = re.sub('_\d+$', '', key)
    renamed_swarms[new_key] = [re.sub('_\d+$', '', x)
                               for x in kept_swarms[key]]

# retrieve the full list of reads from the mothur names dict
names_output = []
for key in renamed_swarms:
    names_list = ''
    for read_name in renamed_swarms[key]:
        if mothur_name_dict[read_name]:
            names_list += mothur_name_dict[read_name] + ','
    names_output.append([key, names_list.rstrip(',')])

# write as names file
with open(output_names, 'w', newline='') as mothur:
    csvwriter = csv.writer(mothur, delimiter=' ')
    csvwriter.writerows(names_output)

# filter fasta by mothur_output
fasta_input = SeqIO.parse(dereplicated_fasta, 'fasta')
records = [x for x in fasta_input]
output_records = [x for x in records if x.id in renamed_swarms.keys()]
SeqIO.write(output_records, output_fasta, 'fasta')
