#!/usr/bin/env Rscript

# usage:
# Rscript /path/to/make_twosample_mle_spp.R sample_bam control_bam chromsizes_file mle_output_file"

## LIBRARIES

require(Rcpp)
require(data.table)
require(spp)
require(rtracklayer)

print("Hello World")

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