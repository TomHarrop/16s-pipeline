#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# hist_files
NameIhistFile <- function(x){
    sn <- strsplit(basename(x), "_")[[1]][1]
    fn <- strsplit(dirname(x), "/")[[1]][c(2, 3)]
    paste(c(fn, sn), collapse = "_")
}

ihist_files <- list.files("output",
           pattern = "merge-ihist.txt",
           full.names = TRUE,
           recursive = TRUE)
names(ihist_files) <- sapply(ihist_files, NameIhistFile)

hist_data_list <- lapply(ihist_files,
                         fread,
                         skip = 6,
                         col.names = c("length", "count"))
hist_data <- rbindlist(hist_data_list, idcol = "samplename")
saveRDS(hist_data, "hist_data.Rds")
ggplot(hist_data, aes(x = length, y = count)) +
    facet_wrap(~samplename, ncol = 3) +
    scale_y_log10() +
    # geom_vline(xintercept = c(200, 218)) +
    geom_col()
