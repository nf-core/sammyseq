#!/usr/bin/env Rscript

# usage:
# Rscript /path/to/make_twosample_mle_spp.R sample_bam control_bam chromsizes_file mle_output_file"

## LIBRARIES

require(Rcpp)
require(data.table)
require(spp)
require(rtracklayer)

################################################
## PARAMS
################################################

# ip_file <- args[1]
# input_file <- args[2]
# chromsizes_file <- args[3]
# mle_output_file <- args[4]

ip_file <- "${bam1}"
input_file <- "${bam2}"
chromsizes_file <- "${chromsizes_file}"
mle_output_file <- "${bam1.baseName}VS${bam2.baseName}_mle.txt"

################################################
# DEFAULT PARAMETERS
################################################

remove_anomalies <- FALSE
debug_mode <- TRUE

################################################
# CORE PROCESSES
################################################

if (debug_mode){
    print(c("ip_file=",ip_file))
    print(c("input_file=",input_file))
    print(c("chromsizes_file=",chromsizes_file))
    print(c("mle_output_file",mle_output_file))
}

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

##sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
##print(sessionInfo())
##sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

# r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
# deseq2.version <- as.character(packageVersion('DESeq2'))

# writeLines(
#     c(
#         '"${task.process}":',
#         paste('    r-base:', r.version),
#         paste('    bioconductor-deseq2:', deseq2.version)
#     ),
# 'versions.yml')

################################################
################################################
################################################
################################################