#!/usr/bin/env Rscript

library(data.table)

# snakemake files
precluster_names <- snakemake@input[["precluster_names"]]
kept_otus_file <- snakemake@output[["kept_otus"]]
count_table_file <- snakemake@output[["count_table"]]
abundance_table_file <- snakemake@output[["abundance_table"]]
filter_file <- snakemake@output[["filter_file"]]

# read the long output
precluster_names <- "output/V6-7/convert_for_gutfilter/long_table.tab"

# generate a sample name from each read
read_names <- fread(precluster_names)
read_names[, read_sample_name := gsub("\\|.+$", "", read_name)]
number_of_samples <- read_names[, length(unique(read_sample_name))]

# define filters
overall_community_thr <- 0.035
sample_5pc_thr <- as.integer(
    round(max(2, 0.05 * number_of_samples), 0))
sample_2pc_thr <- as.integer(
    round(max(2, 0.02 * number_of_samples), 0))
abundance_2pc_thr <- 0.001

# count total reads per sample
total_reads_per_sample <- read_names[, .(
    total_sample_reads = length(unique(read_name))),
    by = read_sample_name]

# count total reads per OTU
total_reads_per_otu <- read_names[, .(
    total_otu_reads = length(unique(read_name))),
    by = otu_id]

# count reads per sample
reads_per_sample <- read_names[, .(
    reads_per_sample = length(unique(read_name))),
    by = .(otu_id, read_sample_name)]

# merge total reads per sample and total reads per otu
reads_per_sample_with_total_per_sample <- merge(
    reads_per_sample, total_reads_per_sample, by = "read_sample_name")
totals <- merge(reads_per_sample_with_total_per_sample,
                total_reads_per_otu, by = "otu_id")

# calculate fractions
totals[, sample_relative_abundance := reads_per_sample / total_sample_reads]

# run filters
totals[, pass_community_thr := 
           sample_relative_abundance >= overall_community_thr]
totals[, pass_2pc_thr := 
           sum(sample_relative_abundance >= abundance_2pc_thr) >=
           sample_2pc_thr,
       by = .(otu_id)]
totals[, pass_5pc_thr := 
           sum(reads_per_sample > 0) > sample_5pc_thr,
       by = otu_id]

# list of kept OTUs
totals[, pass_any_filter :=
           sum(pass_community_thr, pass_2pc_thr, pass_5pc_thr) > 0,
       by = otu_id]
kept_otus <- totals[pass_any_filter == TRUE, unique(otu_id)]

# generate count and abundance tables
filtered_totals <- totals[otu_id %in% kept_otus]
count_table <- dcast(filtered_totals,
      otu_id ~ read_sample_name,
      value.var = "reads_per_sample",
      fill = 0)

abundance_table <- dcast(filtered_totals,
      otu_id ~ read_sample_name,
      value.var = "sample_relative_abundance",
      fill = 0)

# write output
fwrite(totals, filter_file)
fwrite(list(kept_otus), kept_otus_file)
fwrite(count_table, count_table_file)
fwrite(abundance_table, abundance_table_file)
