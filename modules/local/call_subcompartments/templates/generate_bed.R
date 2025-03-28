options("scipen"=999)

# Convert HEX colors to RGB format
rgb_str <- function(hex) {
    paste(as.vector(col2rgb(hex)), collapse = ",")
}

# Function to generate BED and BEDGRAPH files for compartments and eigenvectors
generate_bed_and_bedgraph_files <- function(chromosomes, file_suffix, output_dir, col_names, output_ext) {

    # Ensure the output directory exists
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Read database and extract patient names
    tracks_db <- fread(file = single_dbfile, data.table = FALSE, header = TRUE)
    patients <- unique(tracks_db$Patient_name)

    for (name in patients) {
        file_list <- lapply(chromosomes, function(chr) {
            file_path <- paste0(out_dir, name, "_", chr, file_suffix)
            if (file.exists(file_path)) {
                read.table(file_path, header = TRUE, sep = "\t")
            } else {
                return(NULL)
            }
        })

        file_list <- file_list[!sapply(file_list, is.null)]
        if (length(file_list) == 0) next

        # Combine all chromosome data for the patient
        granges_sample <- do.call("rbind", file_list)
        granges_sample$strand <- "."
        granges_sample$zero <- 0

        # Ensure correct column order with duplicate start and end for coloring
        if ("pc1" %in% col_names) {
            # Eigenvector BEDGRAPH format
            granges_sample <- granges_sample[, c('seqnames', 'start', 'end', 'pc1', 'subcomps_vect', 'zero', 'strand', 'start', 'end', 'subcolor_vect')]
        } else {
            # Compartment BED format
            granges_sample <- granges_sample[, c('seqnames', 'start', 'end', 'subcomps_vect', 'zero', 'strand', 'start', 'end', 'subcolor_vect')]
        }

        # Convert HEX colors to RGB
        hex_cols <- unique(granges_sample$subcolor_vect)
        rgb_conv <- setNames(lapply(hex_cols, rgb_str), hex_cols)

        for (i in hex_cols) {
            granges_sample$subcolor_vect <- gsub(i, rgb_conv[[i]], granges_sample$subcolor_vect)
        }

        # Remove any NA rows to avoid errors in BED/BEDGRAPH files
        granges_sample <- na.omit(granges_sample)

        # Define output file path
        file_name <- paste0(output_dir, name, "_", sub(".tsv", "", file_suffix), output_ext)
        header_bedfile <- paste0('track name="', name, '" description="', name, ' (Emission ordered)" visibility=1 itemRgb="On"')

        # Print file name and row count
        print(paste("Saving:", file_name))

        # Write BED header
        fileConn <- file(file_name)
        writeLines(header_bedfile, fileConn)
        close(fileConn)

        # Write BED or BEDGRAPH data
        write.table(granges_sample, append = TRUE, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    }
}

# Function that generates BED files for compartments and BEDGRAPH files for eigenvectors
generate_compartment_bed_and_bedgraph <- function() {

    # Generate BED files for compartments
    generate_bed_and_bedgraph_files(
        chromosomes = chromosomes,
        file_suffix = "_compartments_merged.tsv",
        output_dir = bed_output_dir,
        col_names = c('seqnames', 'start', 'end', 'subcomps_vect', 'zero', 'strand', 'start', 'end', 'subcolor_vect'),
        output_ext = ".bed"
    )

    # Generate BEDGRAPH files for eigenvectors
    generate_bed_and_bedgraph_files(
        chromosomes = chromosomes,
        file_suffix = "_compartments_eigenvector.tsv",
        output_dir = paste0(work_dir, "bedgraph_eigen_compartment/"),
        col_names = c('seqnames', 'start', 'end', 'pc1', 'subcomps_vect', 'zero', 'strand', 'start', 'end', 'subcolor_vect'),
        output_ext = ".bedgraph"
    )
}
