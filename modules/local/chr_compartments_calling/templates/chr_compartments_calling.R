#!/usr/bin/env Rscript

suppressMessages({
    require(optparse)
    require(data.table)
    require(doParallel)
    require(GenomicRanges)
    require(factoextra)
    require(parallel)
    require(rtracklayer)
    require(ggplot2)
    require(R.utils)
    require(R.oo)
    require(R.methodsS3)
    require(maptools)
    require(patchwork)
    require(sp)
    require(dendextend)
    require(Gviz)
    require(CALDER)
    require(ape)
    require(fitdistrplus)
    require(igraph)
    require(Matrix)
    require(rARPACK)
    require(fields)
    require(strawr)
})

## https://community.seqera.io/t/source-another-r-script-in-the-bin-directory/1059
path <- Sys.getenv("PATH") |> strsplit(":")
bin_path <- tail(path[[1]], n=1)
source(file.path(bin_path, "compartmentalization_analysis_functions.R"))

options(ucscChromosomeNames=FALSE)

## Get parameters from Nextflow template
patient    <- "${patient}"
chromosome <- "${meta.chromosome}"
binsize    <-  as.integer(${binsize})
input_file <- "${csv}"
chrom_bed  <- "${chr_bed}"
gene_gtf   <- "${gtf}"

## Load data
comp_df <- fread(input_file, data.table = FALSE)
comp_df_replica <- comp_df[comp_df\$Patient_name == patient, ]

bins_gr <- import(chrom_bed, format = "BED")
genes_gr <- import(gene_gtf, format = "GTF")

sub2_colors <- c("B" = "#4575b4", "A" = "#d73027")
subs_file <- paste0(patient, "_compartment___", chromosome, '_', binsize, ".Rdata")

# Run SAMMY
sub_objs <- call_subcompartments_sammy(
    patients = patient,
    tracks_db = comp_df_replica,
    bins_gr = bins_gr,
    subs_file = subs_file,
    binsize = binsize,
    chr = chromosome,
    genes_gr = genes_gr,
    keeping_bins1 = "all",
    sublevel = "sub.2",
    sub_colors = sub2_colors
)

generate_files(sub_objs, chromosome)

# Write versions
r_version <- paste(R.Version()[c("major", "minor")], collapse = ".")
writeLines(c(
    '"${task.process}":',
    paste0('    r-base: "', r_version, '"')
), "versions.yml")