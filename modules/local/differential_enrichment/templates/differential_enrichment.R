#!/usr/bin/env Rscript

suppressMessages({
    library(data.table)
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(preprocessCore)
    library(effsize)
    library(dplyr)
    library(GenomicFeatures)
    library(BSDA)
})

source("${projectDir}/bin/differential_solubility_functions.R")

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

opt <- list(
    test_group = "${test_group}",
    ref_group = "${ref_group}",
    contrast_name = "${contrast_name}",
    samplesheet = "${samplesheet}",
    genome_bins = "${genome_bins}",
    binsize = "${binsize}",
    comparison = "${comparison}",
    solubility_threshold = "${solubility_threshold}",
    gtf = "${gtf}"
)

################################################
## PARAMETERS                                 ##
################################################

cat("Processing group comparison:", opt\$test_group, "vs", opt\$ref_group, "\\n")

compare_groups <- opt\$contrast_name
samplesheet_file <- opt\$samplesheet
genome_bins_file <- opt\$genome_bins
binsize_param <- opt\$binsize
comparison_param <- opt\$comparison
threshold <- opt\$solubility_threshold
gtf_file <- opt\$gtf

# Define the summary file variable for sink
summary_file <- paste0(opt\$contrast_name, '_analysis_summary.txt')

# Define threshold for binselector
ths <- as.numeric(threshold)

selected_ratios <- if (grepl(",", comparison_param)) trimws(strsplit(comparison_param, ",")[[1]]) else trimws(comparison_param)

comp_db <- fread(samplesheet_file, data.table = FALSE, header = TRUE)
bins_gr <- import(genome_bins_file, format = "BED")
ratio_col <- which(colnames(comp_db) == 'ratio')
file_col  <- which(colnames(comp_db) == 'file')
id_col    <- which(colnames(comp_db) == 'experimental_id')
group_col <- which(colnames(comp_db) == 'sample_group')

################################################
## FUNCTIONS                                  ##
################################################

# Check if specified groups exist in the data
available_groups <- unique(comp_db[['sample_group']])
specified_groups <- c(opt\$test_group, opt\$ref_group)
missing_groups <- setdiff(specified_groups, available_groups)

if (length(missing_groups) > 0) {
    stop("ERROR: Groups [", paste(missing_groups, collapse = ", "), "] not found in data!\\n",
        "Available groups: [", paste(available_groups, collapse = ", "), "]\\n")
}

# metadata addition function
add_metadata <- function(df, comparison_name, current_ratio, fraction, direction) {
    df[['comparison']] <- comparison_name
    df[['ratio']] <- current_ratio
    df[['fraction']] <- fraction
    df[['direction']] <- direction
    return(df)
}

# save bins data to CSV in a structured way
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
        cat("No", file_suffix, "data found for", comparison_name, "\\n")
    }
}

# Setup gene annotation
if (gtf_file != "" && file.exists(gtf_file)) {
    tryCatch({
        cat("Setting up gene annotations from:", gtf_file, "\\n")
        final_genes <- setup_gene_annotation(gtf_file)
        assign("final_genes", final_genes, envir = .GlobalEnv)
        cat("Gene annotations ready:", length(final_genes), "genes\\n")
    }, error = function(e) {
        cat("WARNING: GTF parsing failed. Skipping gene analysis\\n")
    })
} else {
    cat("No GTF file provided. Gene analysis will be skipped\\n")
}

################################################
## MAIN ANALYSIS                              ##
################################################

sink(summary_file)

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
    gr1_names <- names(gr1@elementMetadata)
    S4Vectors::mcols(gr1) <- preprocessCore::normalize.quantiles(as.matrix(S4Vectors::mcols(gr1)))
    names(gr1@elementMetadata) <- gr1_names
    allmixeddf_grobj <- GenomicRanges::sort(gr1)
    unique_groups <- unique(Sample_groups)


    pr <- matrix(nrow = 2, ncol = 1)
    colnames(pr) <- compare_groups
    rownames(pr) <- c('ref_group', 'test_group')

    pr[1, 1] <- opt\$ref_group
    pr[2, 1] <- opt\$test_group

    assign("pr", pr, envir = .GlobalEnv)
    g1 <- pr[1, 1]; g2 <- pr[2, 1]
    assign(g1, Sample_names[Sample_groups == g1], envir = .GlobalEnv)
    assign(g2, Sample_names[Sample_groups == g2], envir = .GlobalEnv)

    result <- Bins_selector(
        combination = compare_groups,
        allmixeddf_grobj = allmixeddf_grobj,
        fraction1 = fr1,
        fraction2 = fr2,
        ths = ths
    )

    comparison_name <- compare_groups
    all_bins_data <- list()
    selected_bins_data <- list()

    all_bins_result <- result[[paste0(g2, "_allgr_", g1)]]
    if (!is.null(all_bins_result) && length(all_bins_result) > 0) {
        df_all <- as.data.frame(all_bins_result)
        df_all <- add_metadata(df_all, comparison_name, current_ratio, "all", "all")
        all_bins_data[[1]] <- df_all
    } else {
        cat("No all_bins data available for", comparison_name, "\\n")
    }

    bin_categories <- list(
        list(suffix = paste0(g2, "_", fr1, "_up_", g1),   frac = fr1, dir = "up"),
        list(suffix = paste0(g2, "_", fr1, "_down_", g1), frac = fr1, dir = "down"),
        list(suffix = paste0(g2, "_", fr2, "_up_", g1),   frac = fr2, dir = "up"),
        list(suffix = paste0(g2, "_", fr2, "_down_", g1), frac = fr2, dir = "down")
    )

    for (category in bin_categories) {
        gr <- result[[category[['suffix']]]]
        if (!is.null(gr) && length(gr) > 0) {
            df <- as.data.frame(gr)
            df <- add_metadata(df, comparison_name, current_ratio, category[['frac']], category[['dir']])
            selected_bins_data <- c(selected_bins_data, list(df))
        } else {
            cat("No", category[['dir']], "bins found for fraction", category[['frac']], "in", comparison_name, "\\n")
        }
    }

    save_bins_data(all_bins_data, current_ratio, comparison_name, "all_bins_complete", g1, g2)
    save_bins_data(selected_bins_data, current_ratio, comparison_name, "selected_bins_filtered", g1, g2)

    ## Select only relevant GRanges for BED export
    fr1_up_name <- paste0(fr1, "_up")
    fr1_down_name <- paste0(fr1, "_down")
    fr2_up_name <- paste0(fr2, "_up")
    fr2_down_name <- paste0(fr2, "_down")

    selected_bins_only <- list()
    selected_bins_only[[fr1_up_name]] <- result[[paste0(g2, "_", fr1, "_up_", g1)]]
    selected_bins_only[[fr1_down_name]] <- result[[paste0(g2, "_", fr1, "_down_", g1)]]
    selected_bins_only[[fr2_up_name]] <- result[[paste0(g2, "_", fr2, "_up_", g1)]]
    selected_bins_only[[fr2_down_name]] <- result[[paste0(g2, "_", fr2, "_down_", g1)]]

    # Create regions directory
    dir.create("regions", showWarnings = FALSE)

    # Generate BED files for each category
    base_name <- paste0(g2, "vs", g1, "_", current_ratio)
    bed_categories <- c(fr1_up_name, fr1_down_name, fr2_up_name, fr2_down_name)

    for (bed_cat in bed_categories) {
        gr_touse <- selected_bins_only[[bed_cat]]
        if (!is.null(gr_touse) && length(gr_touse) > 0) {
            bed_filename <- paste0("regions/", base_name, "_", bed_cat, "_regions.bed")
            write.table(as.data.frame(gr_touse)[,c(1,2,3)],
                        bed_filename,
                        quote = FALSE,
                        sep = "\\t",
                        row.names = FALSE,
                        col.names = FALSE
                        )
            cat("Saved BED file:", bed_filename, "with", length(gr_touse), "regions\\n")
        } else {
            cat("No regions found for", bed_cat, "in", comparison_name, "\\n")
        }
    }

    # Write gene files if provided
    if ("genes" %in% names(result)) {
        gene_results <- result[["genes"]]
        if (length(gene_results) > 0) {

            dir.create("genes", showWarnings = FALSE)

            for (gene_category in names(gene_results)) {
                gene_list <- gene_results[[gene_category]]
                if (length(gene_list) > 0) {
                    cat("Writing", length(gene_list), "genes for", gene_category, "\\n")

                    output_file <- paste0("genes/", current_ratio, "_", gene_category, "_genes.txt")
                    write.table(gene_list, output_file,
                                quote = FALSE, sep = "\\t",
                                row.names = FALSE, col.names = FALSE)
                }
            }
        } else {
            cat("No gene analysis results available for", comparison_name, "\\n")
        }
    } else {
        cat("Gene analysis was not performed for", comparison_name, "(final_genes not available)\\n")
    }

    # Add quantile normalized bins to result for saving
    result[["quantile_normalized_bins"]] <- allmixeddf_grobj

    # Save complete analysis results as RDS
    dir.create("rdata", showWarnings = FALSE)
    output_complete_rds <- paste0("rdata/", current_ratio, "_", compare_groups, "_analysis.rds")
    saveRDS(result, file = output_complete_rds)
    cat("Saved complete analysis R object with normalized data:", output_complete_rds, "\\n")
}

sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
writeLines(
    c(
    '"${task.process}":',
    paste('    r-base:', r.version),
    paste('    r-data.table:', as.character(packageVersion('data.table'))),
    paste('    bioconductor-rtracklayer:', as.character(packageVersion('rtracklayer'))),
    paste('    bioconductor-genomicranges:', as.character(packageVersion('GenomicRanges'))),
    paste('    bioconductor-genomeinfodb:', as.character(packageVersion('GenomeInfoDb'))),
    paste('    bioconductor-preprocesscore:', as.character(packageVersion('preprocessCore'))),
    paste('    r-effsize:', as.character(packageVersion('effsize'))),
    paste('    r-dplyr:', as.character(packageVersion('dplyr'))),
    paste('    bioconductor-genomicfeatures:', as.character(packageVersion('GenomicFeatures'))),
    paste('    r-bsda:', as.character(packageVersion('BSDA')))
    ),
'versions.yml')
