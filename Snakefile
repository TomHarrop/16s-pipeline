#!/usr/bin/env python3

import os
import pandas


#############
# FUNCTIONS #
#############

def trim_input_resolver(wildcards):
    # get all the matching filenames
    filenames = sorted(set(
        sample_df[(sample_df['Mean Quality Score (PF)'] > 30) &
                  (sample_df.amplicon == wildcards.amplicon) &
                  (sample_df.sampletype == wildcards.sampletype) &
                  (sample_df.samplename == wildcards.samplename)].filename))
    sample_files = []
    for fn in filenames:
        for fq in fastq_files:
            if fn in fq:
                sample_files.append(fq)
    # separate into R1 and R2
    input_files = {}
    input_files['R1'] = sorted(set(x for x in sample_files if 'R1' in x))
    input_files['R2'] = sorted(set(x for x in sample_files if 'R2' in x))

    return input_files


###########
# GLOBALS #
###########

read_dir = 'data/NZGL01875'
sample_key = 'data/nzgl01875_seq_summary.csv'

#########
# SETUP #
#########

# mung the sample key
sample_df = pandas.read_csv(sample_key)
sample_df['samplename'] = sample_df['Sample Ref'].str.split("-").str.get(0)
sample_df['sampletype'] = sample_df['Sample Ref'].str.split("-").str.get(1)
sample_df['amplicon'] = sample_df['Library Type'].str.replace("/", "-")
sample_df.rename(
    columns={
        'OGBF ID': 'filename'
    }, inplace=True)

# get a list of fastq files
read_dir_files = [(dirpath, filenames) for (dirpath, dirnames, filenames)
                  in os.walk(read_dir, followlinks=True)]
fastq_files = []
for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            fastq_files.append(os.path.join(dirpath, filename))

# generate the target paths
trim_target_paths = sorted(
    set('{0}/{1}/{2}'.format(
            x[1].amplicon,
            x[1].sampletype,
            x[1].samplename)
        for x in sample_df.iterrows()))
all_samples = sorted(set(sample_df.samplename))


#########
# RULES #
#########
rule all:
    input:
        expand(('output/{amplicon}/{sampletype}/'
                'convert_for_gutfilter/precluster.fasta'),
               amplicon=['V1-2', 'V6-7'],
               sampletype=['DNA', 'RNA'])

rule trim_merge:
    input:
        unpack(trim_input_resolver)
    output:
        ('output/{amplicon}/{sampletype}/trim_merge/'
         '{samplename}.fastq.gz')
    params:
        stats = ('output/{amplicon}/{sampletype}/trim_merge/'
                 '{samplename}_trim-stats.txt'),
        ihist = ('output/{amplicon}/{sampletype}/trim_merge/'
                 '{samplename}_merge-ihist.txt'),
        adaptors = ('output/{amplicon}/{sampletype}/trim_merge/'
                    '{samplename}_merge-adaptors.txt'),
        ref = 'data/primers.fasta'
    log:
        bbduk = ('output/{amplicon}/{sampletype}/trim_merge/'
                 '{samplename}_bbduk.log'),
        bbmerge = ('output/{amplicon}/{sampletype}/trim_merge/'
                   '{samplename}_bbmerge.log')
    threads:
        8
    shell:
        'bin/bbmap/bbduk.sh '
        'in={input.R1} '
        'in2={input.R2} '
        'out=stdout.fastq '
        'stats={params.stats} '
        'forcetrimright=199 '
        'k=27 '
        'mink=10 '
        'ktrim=l tpe tbo '
        'editdistance=2 '
        'editdistance2=2 '
        'ref={params.ref} '
        '2>> {log.bbduk} '
        ' | '
        'bin/bbmap/bbmerge.sh '
        'in=stdin.fastq '
        'out=stdout.fastq '
        'outa={params.adaptors} '
        'maxlength=219 '
        'mininsert=200 '
        'adapter={params.ref} '
        'ihist={params.ihist} '
        'strict=t '
        '2> {log.bbmerge}'
        ' | '
        'bin/bbmap/bbduk.sh '
        'in=stdin.fastq '
        'maxns=0 '
        'literal='
        'AAAAAAAAAA,'
        'CCCCCCCCCC,'
        'GGGGGGGGGG,'
        'TTTTTTTTTT '
        'maskmiddle=f '
        'out={output} '
        'ziplevel=9 '
        'forcetrimright=199 '
        '2>> {log.bbduk} '

rule rename_reads:
    input:
        ('output/{amplicon}/{sampletype}/trim_merge/'
         '{samplename}.fastq.gz')
    output:
        ('output/{amplicon}/{sampletype}/rename_reads/'
         '{samplename}.fasta')
    params:
        samplename = '{samplename}'
    script:
        'src/rename_fastq_header.py'

rule merge_samples:
    input:
        expand('output/{{amplicon}}/{{sampletype}}/rename_reads/'
               '{sample}.fasta',
               sample=all_samples)
    output:
        'output/{amplicon}/{sampletype}/merge_samples/all.fasta'
    threads:
        1
    shell:
        'cat {input} > {output}'

# Dereplicate using mothur
rule dereplicate:
    input:
        'output/{amplicon}/{sampletype}/merge_samples/all.fasta'
    output:
        fa = temp('output/{amplicon}/{sampletype}/dereplicate/all.fasta'),
        derep = 'output/{amplicon}/{sampletype}/dereplicate/all.unique.fasta',
        names = 'output/{amplicon}/{sampletype}/dereplicate/all.names'
    threads:
        1
    params:
        wd = 'output/{amplicon}/{sampletype}/dereplicate/'
    log:
        'output/{amplicon}/{sampletype}/dereplicate/mothur.log'
    shell:
        'cp {input} {output.fa} ; '
        'bash -c \''
        'cd {params.wd} || exit 1 ; '
        'mothur "'
        '#unique.seqs(fasta=all.fasta)" '
        '\' &> {log}'

# Reformat mothur’s output for swarm:
rule reformat_for_swarm:
    input:
        'output/{amplicon}/{sampletype}/dereplicate/all.unique.fasta'
    output:
        'output/{amplicon}/{sampletype}/reformat_for_swarm/all.unique.fasta'
    params:
        names = 'output/{amplicon}/{sampletype}/dereplicate/all.names'
    threads:
        1
    script:
        'src/reheader_dereplicated_fasta.py'

# Cluster using swarm:
rule swarm:
    input:
        'output/{amplicon}/{sampletype}/reformat_for_swarm/all.unique.fasta'
    output:
        'output/{amplicon}/{sampletype}/swarm/all.unique.swarm'
    threads:
        8
    log:
        'output/{amplicon}/{sampletype}/swarm/swarm.log'
    shell:
        'swarm '
        '-d 6 '
        '-t {threads} '
        '-l {log} '
        '-o {output} '
        '{input}'

# Convert the output for gutfilter
rule convert_for_gutfilter:
    input:
        mothur_names = 'output/{amplicon}/{sampletype}/dereplicate/all.names',
        swarm_results = 'output/{amplicon}/{sampletype}/swarm/all.unique.swarm',
        dereplicated_fasta = ('output/{amplicon}/{sampletype}/dereplicate/'
                              'all.unique.fasta')
    output:
        names = ('output/{amplicon}/{sampletype}/'
                 'convert_for_gutfilter/precluster.names'),
        fasta = ('output/{amplicon}/{sampletype}/'
                 'convert_for_gutfilter/precluster.fasta')
    threads:
        1
    script:
        'src/reformat_swarm_output.py'
