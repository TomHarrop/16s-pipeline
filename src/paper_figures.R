#!/usr/bin/env Rscript

library(DESeq2)
library(data.table)
library(phyloseq)
library(ggplot2)
library(scales)

#############
# FUNCTIONS #
#############

SplitTaxonomy <- function(x) {
    no_nums <- gsub("\\([[:digit:]]+\\)", "", x)
    strsplit(no_nums, ";")
}


###########
# GLOBALS #
###########

outdir <- "paper_output"

# files
count_file <- "output/V6-7/gutfilter/count_table.txt"
taxonomy_file <- "output/V6-7/annotate_otus/keptotus.seed_v128.wang.taxonomy"
tree_file <- "output/V6-7/tree/keptotus.phylip.taxonomy.tre"

tax_levels <- c("Domain",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus")

samplename_grep <- "^X?([[:alnum:]]+)_DNA"

#########
# SETUP #
#########

if(!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
}

# generate sample data
sample_table <- data.table(samplename = colnames(count_matrix))
sample_table[, c("individual", "sampletype") := tstrsplit(samplename, "_")]
sd <- sample_data(data.frame(
    sample_table[sampletype == "DNA", .(samplename, individual)],
    row.names = "samplename"))
sd$samplename <- rownames(sd)

# generate OTU table
count_table <- fread(count_file)
keep_cols <- grep("DNA$", names(count_table), value = TRUE)
filtered_otus <- melt(count_table[, c("otu_id", keep_cols),
                                  with = FALSE],
                      id = "otu_id")[,
                                     sum(value) > 0,
                                     by = .(otu_id, variable)][
                                         V1 == TRUE, unique(otu_id)]
filtered_count_table <- count_table[otu_id %in% filtered_otus,
                                    c("otu_id", keep_cols),
                                    with = FALSE]
count_matrix <- as.matrix(data.frame(filtered_count_table,
                                     row.names = "otu_id"))
otu <- otu_table(count_matrix, taxa_are_rows = TRUE)

# generate taxonomy matrix
taxonomy_table <- fread(taxonomy_file,
                        header = FALSE,
                        col.names = c("otu_id", "taxonomy"))
taxonomy_by_otu <- taxonomy_table[, transpose(SplitTaxonomy(taxonomy)),
                                  by = otu_id]
setnames(taxonomy_by_otu, paste0("V", c(1:6)), tax_levels)
tax_matrix <- as.matrix(data.frame(taxonomy_by_otu, row.names = "otu_id"))
tax <- tax_table(tax_matrix)

# generate tree
tree <- phyloseq::read_tree(tree_file)

# generate the phyloseq object
physeq <- phyloseq(otu, tax, sd, tree)

# convert to DESeq2
dds <- phyloseq_to_deseq2(physeq, ~ 1)
dds <- estimateSizeFactors(dds)

# transform read counts
vst <- varianceStabilizingTransformation(dds, blind = FALSE)

#################
# CHAO RICHNESS #
#################

chao <- data.table(
    estimate_richness(physeq, measures = "Chao1"),
    keep.rownames = TRUE
)
setnames(chao, "rn", "Sample")
chao[, Sample := gsub(samplename_grep, "\\1", Sample)]
fwrite(chao, paste0(outdir, "/alpha_diversity.csv"))

#########################
# OTU BAR PLOT from VST #
#########################

vst_counts_wide <- data.table(assay(vst), keep.rownames = TRUE)
vst_counts <- melt(vst_counts_wide,
                   id = 'rn',
                   variable.name = "Sample",
                   value.name = "normalized_counts")

vst_counts[, Sample := gsub(samplename_grep, "\\1", Sample)]
setnames(vst_counts, "rn", "OTU")

vst_counts_tax <- merge(vst_counts,
                        taxonomy_by_otu,
                        by.x = "OTU",
                        by.y = "otu_id",
                        all.x = TRUE)

vst_pd <- vst_counts_tax[!grepl("unclassified", Genus),
                         sum(normalized_counts),
                         by = .(Sample, Phylum, Genus)]

ggplot(vst_pd, aes(x = Genus, y = V1, fill = Phylum)) +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    facet_grid(Sample ~ Phylum, scales = "free_x", space = "free_x") +
    geom_col(position = "dodge", width = 0.5, size = 0.5)


############################
# OTU BAR PLOT from COUNTS #
############################

counts_wide <- counts(dds, normalized = TRUE)
norm_counts <- melt(data.table(counts_wide, keep.rownames = TRUE),
                    id = 'rn',
                    variable.name = "Sample",
                    value.name = "normalized_counts")
norm_counts[, Sample := gsub(samplename_grep, "\\1", Sample)]
setnames(norm_counts, "rn", "OTU")

norm_counts_tax <- merge(norm_counts,
                         taxonomy_by_otu,
                         by.x = "OTU",
                         by.y = "otu_id",
                         all.x = TRUE)

norm_counts_by_genus <- norm_counts_tax[, .(total_reads = sum(normalized_counts,
                                                              na.rm = TRUE)),
                                        by = .(Sample, Genus, Phylum)]
norm_counts_by_genus[, sum_by_genus := sum(total_reads), by = Genus]
norm_counts_by_genus[, sum_by_sample := sum(total_reads), by = Sample]
norm_counts_by_genus[, pct_reads := total_reads * 100 / sum_by_sample]
setorder(norm_counts_by_genus, Phylum, -sum_by_genus)
norm_counts_by_genus[, Genus := factor(Genus, levels = unique(Genus))]

plot_genus <- norm_counts_by_genus[!grepl("unclassified", Genus),
                                   any(total_reads > 10),
                                   by = Genus][
                                       V1 == TRUE, unique(Genus)]
plot_otu <- norm_counts_tax[Genus %in% as.character(plot_genus), unique(OTU)]

bp_pd <- norm_counts_by_genus[Genus %in% plot_genus]

gp <- ggplot(bp_pd[total_reads > 0],
             aes(x = Genus, y = total_reads, fill = Phylum)) +
    theme_minimal(base_size = 12) +
    scale_fill_brewer(palette = "Set1") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    ylab("Abundance") +
    facet_grid(Sample ~ Phylum, scales = "free_x", space = "free_x") +
    geom_col(position = "dodge", width = 0.5, size = 0.5)

ggsave(paste0(outdir, "/otu_abundance.pdf"),
       gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in",
       family = "Helvetica")

#################
# BRAY DISTANCE #
#################

# calculate sample order
sn_dist <- distance(physeq, method = "bray")
sn_hc <- hclust(sn_dist, method = "average")
plot(sn_hc)
x_sample_order <- sn_hc$labels[sn_hc$order]
sample_order <- gsub("^X", "", x_sample_order)

# distance table
dist_labels <- attr(sn_dist, "Labels")
attr(sn_dist, "Labels") <- gsub(samplename_grep, "\\1", dist_labels)
bray <- data.table(as.matrix(sn_dist), keep.rownames = TRUE)
setnames(bray, "rn", "")
fwrite(bray, paste0(outdir, "/bray_sample_distance.csv"))

###############
# OTU HEATMAP #
###############

# scale the VST
scaled_counts <- t(scale(t(assay(vst))))[plot_otu, ]

# otu order for scaled counts
hc_scaled <- hclust(dist(scaled_counts,
                         method = "minkowski"))
otu_order_scaled <- hc_scaled$labels[hc_scaled$order]
hc_cuts <- cutree(hc_scaled, k = 4)
hc_groups <- data.table(otu_id = names(hc_cuts),
                        group = as.integer(hc_cuts))

# set up plot data
sc_dt_wide <- data.table(data.frame(scaled_counts),
                         keep.rownames = TRUE)
sc_long <- melt(sc_dt_wide,
                id.vars = "rn",
                variable.name = "sample_name",
                value.name = "scaled_counts")
setnames(sc_long, "rn", "otu_id")

# add taxonomy names
sc_pd <- merge(merge(sc_long, taxonomy_by_otu),
               hc_groups)[Genus %in% plot_genus]

# sample name order
sc_pd[, sample_name := factor(gsub(samplename_grep,
                                   "\\1",
                                   sample_name),
                              levels = gsub(samplename_grep,
                                            "\\1",
                                            x_sample_order))]

# OTU labels
sc_pd[, otu_id := factor(otu_id, levels = otu_order_scaled)]
setorder(sc_pd, otu_id)
y_labs <- unique(sc_pd,
                 by = c("otu_id",
                        "Genus"))[, paste0(otu_id, " (", Genus, ")")]
names(y_labs) <- unique(sc_pd, by = c("otu_id", "Genus"))[, otu_id]

# OTU label colour
sc_pd[, Phylum := factor(Phylum, levels = sort(unique(Phylum)))]
y_colour <- RColorBrewer::brewer.pal(5, "Set1")
names(y_colour) <- sc_pd[, levels(Phylum)]

sc_pd[, label_colour := y_colour[as.character(Phylum)], by = .(otu_id)]
col_order <- unique(sc_pd, by = c("otu_id"))$label_colour

# plot
gp <- ggplot(sc_pd, aes(x = sample_name,
                        y = otu_id,
                        fill = scaled_counts,
                        colour = Phylum)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(colour = col_order)) +
    xlab(NULL) + ylab("OTU (Genus)") +
    scale_fill_gradientn(colours = hs,
                         guide = guide_colourbar(title = "Scaled abundance")) +
    scale_y_discrete(labels = y_labs,
                     expand = c(0,0)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_colour_manual(values = y_colour) +
    geom_text(size = 0, label = "") +
    geom_raster() +
    guides(colour = guide_legend(
        title = "Phylum",
        override.aes = list(size = 5)))

ggsave(paste0(outdir, "/otu_heatmap.pdf"),
       gp,
       device = cairo_pdf,
       width = 10,
       height = 7.5,
       units = "in",
       family = "Helvetica")

