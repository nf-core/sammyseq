#!/usr/bin/env Rscript

library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)
library(preprocessCore)
library(effsize)        
library(dplyr)
library(GenomicFeatures)

source("${projectDir}/bin/differential_solubility_functions.R")

samplesheet_file <- "${samplesheet}"
genome_bins_file <- "${genome_bins}"
binsize_param    <- "${binsize}"
comparison_param <- "${comparison}"
compare_groups <- "${compare_groups}"
threshold <- "${solubility_threshold}"
gtf_file <- "${gtf}"

# Define name for the text report
summary_file <- 'analysis_summary.txt'

# Define threshold for binselector
ths <- as.numeric(threshold)

selected_ratios <- if (grepl(",", comparison_param)) trimws(strsplit(comparison_param, ",")[[1]]) else trimws(comparison_param)

comp_db <- fread(samplesheet_file, data.table = FALSE, header = TRUE)
bins_gr <- import(genome_bins_file, format = "BED")
ratio_col <- which(colnames(comp_db) == 'ratio')
file_col  <- which(colnames(comp_db) == 'file')
id_col    <- which(colnames(comp_db) == 'experimental_id')
group_col <- which(colnames(comp_db) == 'sample_group')

# Function to add metadata columns to a dataframe
add_metadata <- function(df, comparison_name, current_ratio, fraction, direction) {
    df[['comparison']] <- comparison_name
    df[['ratio']] <- current_ratio
    df[['fraction']] <- fraction
    df[['direction']] <- direction
    return(df)
}

# Function to save bins data into one CSV with a readable column order
save_bins_data <- function(data_list, current_ratio, comparison_name, file_suffix, g1 = NULL, g2 = NULL) {
    if (length(data_list)) {
        df <- do.call(rbind, data_list)
        if (grepl("all_bins", file_suffix, ignore.case = TRUE)) {
            cols_to_remove <- c('fraction', 'direction')
            available_cols <- colnames(df)
            cols_to_keep <- setdiff(available_cols, cols_to_remove)
            df <- df[, cols_to_keep]
        } else {
            base_cols <- c('seqnames', 'start', 'end', 'ratio', 'comparison', 'fraction', 'direction')
            essential_stats_cols <- c()
            if (!is.null(g1) && !is.null(g2)) {
                essential_stats_cols <- c(
                    paste0(g1, "_serrx2_lower"), paste0(g1, "_serrx2_upper"),
                    paste0(g1, "_mean"), paste0(g1, "_serrX2"),
                    paste0(g2, "_serrx2_lower"), paste0(g2, "_serrx2_upper"),
                    paste0(g2, "_mean"), paste0(g2, "_serrX2"),
                    "delta",
                    "cohen.estimate",
                    "cohen.magnitude"
                    )}
            essential_cols <- c(base_cols, essential_stats_cols)
            available_cols <- colnames(df)
            cols_to_keep <- intersect(essential_cols, available_cols)
            df <- df[, cols_to_keep]
        }
        output_file <- paste0(current_ratio, "_", comparison_name, "_", file_suffix, ".csv")
        write.csv(df, file = output_file, quote = FALSE, row.names = FALSE)
    } else {
        cat("No", file_suffix, "data found for", comparison_name, "\n")
    }
}

sink(summary_file)

# Setup gene annotations if GTF is provided
if (gtf_file != "" && file.exists(gtf_file)) {
    tryCatch({
        cat("Setting up gene annotations from:", gtf_file, "\n")
        final_genes <- setup_gene_annotation(gtf_file, NULL)
        assign("final_genes", final_genes, envir = .GlobalEnv)
        cat("Gene annotations ready:", length(final_genes), "genes\n")
    }, error = function(e) {
        cat("WARNING: GTF parsing failed. Skipping gene analysis\n")
    })
} else {
    cat("No GTF file provided. Gene analysis will be skipped\n")
}

for (current_ratio in selected_ratios) {
    ratio_data    <- comp_db[comp_db[, ratio_col] == current_ratio, ]
    Sample_names  <- ratio_data[, id_col]
    Sample_groups <- ratio_data[, group_col]
    Sample_files  <- ratio_data[, file_col]
    fr_parts <- strsplit(current_ratio, "vs", perl = TRUE)[[1]]
    fr1 <- trimws(fr_parts[1]); fr2 <- trimws(fr_parts[2])

    bws <- import_and_rebin__bw(
        files    = Sample_files,
        bin_list = bins_gr,
        names    = Sample_names,
        genome   = NULL
    )

    bindf <- as.data.frame(bins_gr)[c(1,2,3,4,5)]
    for (sample_name in names(bws)) {
        dftomerge <- as.data.frame(bws[[sample_name]])
        colnames(dftomerge)[6] <- sample_name
        bindf <- merge(bindf, dftomerge[, 1:6], by = c(1,2,3,4,5), sort = FALSE)
    }

    gr1 <- GenomicRanges::makeGRangesFromDataFrame(bindf, keep.extra.columns = TRUE)
    # gr2 <- GenomicRanges::makeGRangesFromDataFrame(bindf, keep.extra.columns = TRUE)
    gr1_names <- names(gr1@elementMetadata)
    S4Vectors::mcols(gr1) <- preprocessCore::normalize.quantiles(as.matrix(S4Vectors::mcols(gr1)))
    # names(gr1@elementMetadata) <- names(gr2@elementMetadata)
    names(gr1@elementMetadata) <- gr1_names
    allmixeddf_grobj <- GenomicRanges::sort(gr1)
    unique_groups <- unique(Sample_groups)

    # Process custom comparisons
    comps <- strsplit(compare_groups, ",")[[1]]
    pr <- matrix(nrow = 2, ncol = length(comps))
    colnames(pr) <- comps
    rownames(pr) <- c('ref_group', 'test_group')
    for (i in seq_along(comps)) {
        parts <- strsplit(comps[i], "vs")[[1]]
        test_group <- trimws(parts[1])
        ref_group <- trimws(parts[2])
        if (!test_group %in% unique_groups) {
            stop("Group '", test_group, "' not found in ", current_ratio, ". Available groups: ", paste(unique_groups, collapse = ", "))
        }
        if (!ref_group %in% unique_groups) {
            stop("Group '", ref_group, "' not found in ", current_ratio, ". Available groups: ", paste(unique_groups, collapse = ", "))
        }
        pr[1, i] <- ref_group
        pr[2, i] <- test_group
    }
    cat("Custom comparisons:", paste(apply(pr, 2, function(x) paste0(x[2], "_vs_", x[1])), collapse = ", "), "\n")

    assign("pr", pr, envir = .GlobalEnv)
    for (i in 1:ncol(pr)) {
        g1 <- pr[1, i]; g2 <- pr[2, i]
        assign(g1, Sample_names[Sample_groups == g1], envir = .GlobalEnv)
        assign(g2, Sample_names[Sample_groups == g2], envir = .GlobalEnv)
    }

    list_groups <- vector("list", ncol(pr))
    names(list_groups) <- comps
    for (i in 1:ncol(pr)) {
        list_groups[[i]] <- Bins_selector(
            combination = comps[i],
            allmixeddf_grobj = allmixeddf_grobj,
            fraction1 = fr1,
            fraction2 = fr2,
            ths = ths
        )
    }
    # names(list_groups) <- apply(pr, 2, function(x) paste(x[1], "vs", x[2], sep = "_"))

    for (i in seq_along(list_groups)) {
        g1 <- pr[1, i]; g2 <- pr[2, i]
        res <- list_groups[[i]]
        comparison_name <- names(list_groups)[i]
        all_bins_data <- list()
        selected_bins_data <- list()

        all_bins_result <- res[[paste0(g2, "_allgr_", g1)]]
        if (!is.null(all_bins_result) && length(all_bins_result) > 0) {
            df_all <- as.data.frame(all_bins_result)
            df_all <- add_metadata(df_all, comparison_name, current_ratio, "all", "all")
            all_bins_data[[1]] <- df_all
        } else {
            cat("No all_bins data available for", comparison_name, "\n")
        }

        bin_categories <- list(
            list(suffix = paste0(g2, "_", fr1, "_up_", g1),   frac = fr1, dir = "up"),
            list(suffix = paste0(g2, "_", fr1, "_down_", g1), frac = fr1, dir = "down"),
            list(suffix = paste0(g2, "_", fr2, "_up_", g1),   frac = fr2, dir = "up"),
            list(suffix = paste0(g2, "_", fr2, "_down_", g1), frac = fr2, dir = "down")
        )

        for (category in bin_categories) {
            gr <- res[[category[['suffix']]]]
            if (!is.null(gr) && length(gr) > 0) {
                df <- as.data.frame(gr)
                df <- add_metadata(df, comparison_name, current_ratio, category[['frac']], category[['dir']])
                selected_bins_data <- c(selected_bins_data, list(df))
            } else {
                cat("No", category[['dir']], "bins found for fraction", category[['frac']], "in", comparison_name, "\n")
            }
        }

        save_bins_data(all_bins_data, current_ratio, comparison_name, "all_bins_complete", g1, g2)
        save_bins_data(selected_bins_data, current_ratio, comparison_name, "selected_bins_filtered", g1, g2)

        # Write gene files if available
        if ("genes" %in% names(res)) {
            gene_results <- res[["genes"]]
            if (length(gene_results) > 0) {
                for (gene_category in names(gene_results)) {
                    gene_list <- gene_results[[gene_category]]
                    if (length(gene_list) > 0) {
                        cat("Writing", length(gene_list), "genes for", gene_category, "\n")
                        
                        output_file <- paste0(current_ratio, "_", gene_category, "_genes.txt")
                        write.table(gene_list, output_file,
                                    quote = FALSE, sep = "\t", 
                                    row.names = FALSE, col.names = FALSE)
                    }
                }
            } else {
                cat("No gene analysis results available for", comparison_name, "\n")
            }
        } else {
            cat("Gene analysis was not performed for", comparison_name, "(final_genes not available)\n")
        }
    }

    # Save complete analysis results as R object
    comparison_names <- gsub(",", "_", compare_groups)
    output_complete_rds <- paste0(current_ratio, "_", comparison_names, "_complete_analysis.rds")
    saveRDS(list_groups, file = output_complete_rds)
    cat("Saved complete analysis R object:", output_complete_rds, "\n")
}

sink()