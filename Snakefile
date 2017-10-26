#!/usr/bin/env python3

import os
import pandas
import pathlib
import shutil

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


def get_full_path(binary):
    which = shutil.which(binary)
    # check if the binary exists
    if not which:
        raise EnvironmentError(
            'Dependency {0} not found in $PATH'.format(binary))
    # get the full path to binary
    binary_path = pathlib.Path(which).resolve()
    return str(binary_path)


###########
# GLOBALS #
###########

read_dir = 'data/NZGL01875'
sample_key = 'data/nzgl01875_seq_summary.csv'

# binaries
mothur = get_full_path('bin/mothur/mothur')
swarm = get_full_path('bin/swarm')

# datafiles
silva_align = pathlib.Path('data/silva/silva.seed_v128.align').resolve()
silva_tax = pathlib.Path('data/silva/silva.seed_v128.tax').resolve()

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
all_samples = sorted(set(sample_df.samplename))

#########
# RULES #
#########
rule all:
    input:
        expand(('output/{amplicon}/tree/'
                'keptotus.phylip.tre'),
               amplicon=['V1-2', 'V6-7'])

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
        samplename = '{samplename}',
        sampletype = '{sampletype}'
    script:
        'src/rename_fastq_header.py'

rule merge_samples:
    input:
        expand('output/{{amplicon}}/{sampletype}/rename_reads/'
               '{samplename}.fasta',
               samplename=all_samples,
               sampletype=['DNA', 'RNA'])
    output:
        'output/{amplicon}/merge_samples/all.fasta'
    threads:
        1
    shell:
        'cat {input} > {output}'

# Dereplicate using mothur
rule dereplicate:
    input:
        'output/{amplicon}/merge_samples/all.fasta'
    output:
        fa = temp('output/{amplicon}/dereplicate/all.fasta'),
        derep = 'output/{amplicon}/dereplicate/all.unique.fasta',
        names = 'output/{amplicon}/dereplicate/all.names'
    threads:
        1
    params:
        wd = 'output/{amplicon}/dereplicate/'
    log:
        'output/{amplicon}/dereplicate/mothur.log'
    shell:
        'cp {input} {output.fa} ; '
        'bash -c \''
        'cd {params.wd} || exit 1 ; '
        '{mothur} "'
        '#unique.seqs(fasta=all.fasta)" '
        '\' &> {log}'

# Reformat mothur’s output for swarm:
rule reformat_for_swarm:
    input:
        'output/{amplicon}/dereplicate/all.unique.fasta'
    output:
        'output/{amplicon}/reformat_for_swarm/all.unique.fasta'
    params:
        names = 'output/{amplicon}/dereplicate/all.names'
    threads:
        1
    script:
        'src/reheader_dereplicated_fasta.py'

# Cluster using swarm:
rule swarm:
    input:
        'output/{amplicon}/reformat_for_swarm/all.unique.fasta'
    output:
        'output/{amplicon}/swarm/all.unique.swarm'
    threads:
        8
    log:
        'output/{amplicon}/swarm/swarm.log'
    shell:
        '{swarm} '
        '-d 6 '
        '-t {threads} '
        '-l {log} '
        '-o {output} '
        '{input}'

# Convert the output for gutfilter
rule convert_for_gutfilter:
    input:
        mothur_names = 'output/{amplicon}/dereplicate/all.names',
        swarm_results = 'output/{amplicon}/swarm/all.unique.swarm',
        dereplicated_fasta = ('output/{amplicon}/dereplicate/'
                              'all.unique.fasta')
    output:
        names = ('output/{amplicon}/'
                 'convert_for_gutfilter/precluster.names'),
        fasta = ('output/{amplicon}/'
                 'convert_for_gutfilter/precluster.fasta'),
        long_table = ('output/{amplicon}/'
                      'convert_for_gutfilter/long_table.tab'),
    threads:
        1
    script:
        'src/reformat_swarm_output.py'

rule gutfilter:
    input:
        precluster_names = ('output/{amplicon}/'
                            'convert_for_gutfilter/long_table.tab')
    output:
        kept_otus = 'output/{amplicon}/gutfilter/kept_otus.txt',
        count_table = 'output/{amplicon}/gutfilter/count_table.txt',
        abundance_table = 'output/{amplicon}/gutfilter/abundance_table.txt',
        filter_file = 'output/{amplicon}/gutfilter/filter_table.txt'
    threads:
        1
    script:
        'src/gut_filter.R'

rule gutfilter_reads:
    input:
        kept_otus = 'output/{amplicon}/gutfilter/kept_otus.txt',
        fasta = ('output/{amplicon}/'
                 'convert_for_gutfilter/precluster.fasta')
    output:
        'output/{amplicon}/gutfilter/keptotus.fasta'
    log:
        'output/{amplicon}/gutfilter/filterbyname.log'
    threads:
        1
    shell:
        'bin/bbmap/filterbyname.sh '
        'in={input.fasta} '
        'out={output} '
        'names={input.kept_otus} '
        'include=t '
        '2> {log}'

# Annotate the OTUs in the table using mothur
rule annotate_otus:
    input:
        fasta = 'output/{amplicon}/gutfilter/keptotus.fasta'
    output:
        fasta = temp('output/{amplicon}/annotate_otus/keptotus.fasta'),
        tax = 'output/{amplicon}/annotate_otus/keptotus.seed_v128.wang.taxonomy'
    params:
        wd = 'output/{amplicon}/annotate_otus'
    log:
        'output/{amplicon}/annotate_otus/mothur.log'
    threads:
        16
    shell:
        'cp {input.fasta} {output.fasta} ; '
        'bash -c \''
        'cd {params.wd} || exit 1 ; '
        '{mothur} "'
        '#classify.seqs(fasta=keptotus.fasta, '
        'template={silva_align}, '
        'taxonomy={silva_tax}, '
        'processors={threads})" '
        '\' &> {log}'

# Make a phylogenetic tree with mothur
rule mothur_tree:
    input:
        fasta = 'output/{amplicon}/gutfilter/keptotus.fasta',
        tax = ('output/{amplicon}/annotate_otus/'
               'keptotus.seed_v128.wang.taxonomy')
    output:
        fasta = temp('output/{amplicon}/tree/keptotus.fasta'),
#        align = 'output/{amplicon}/tree/keptotus.align',
        tax = temp(('output/{amplicon}/tree/'
                    'keptotus.seed_v128.wang.taxonomy')),
        tree = 'output/{amplicon}/tree/keptotus.phylip.tre'
    params:
        wd = 'output/{amplicon}/tree'
    threads:
        16
    log:
        'output/{amplicon}/tree/mothur.log'
    shell:
        'cp {input.fasta} {output.fasta} ; '
        'cp {input.tax} {output.tax} ; '
        'bash -c \''
        'cd {params.wd} || exit 1 ; '
        '{mothur} "'
        '#align.seqs(candidate=keptotus.fasta, '
        'template={silva_align}, '
        'processors={threads}, '
        'flip=t)" ; '
        '{mothur} "#dist.seqs(fasta=keptotus.align, output=lt)" ; '
        '{mothur} "#tree.shared(phylip=keptotus.phylip.dist)" ; '
        '{mothur} "#clearcut(phylip=keptotus.phylip.dist)" ; '
        '{mothur} "#classify.tree(taxonomy=keptotus.seed_v128.wang.taxonomy, '
        'tree=keptotus.phylip.tre)" '
        '\' &> {log} '
