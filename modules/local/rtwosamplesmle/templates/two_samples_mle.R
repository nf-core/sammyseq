#!/usr/bin/env Rscript

# usage:
# Rscript /path/to/make_twosample_mle_spp.R sample_bam control_bam chromsizes_file mle_output_file"

## LIBRARIES

require(Rcpp)
require(data.table)
require(spp)
require(rtracklayer)

################################################
# FUNCTIONS
################################################

# sortbychr <- function(x, chrcol="chr", stcol="start", endcol=NULL, chrorder=paste("chr", c(seq(22), "X", "Y"), sep="")) {
#         if (!(is.null(endcol))) {
#                 x <- x[order(x[,endcol]),]
#         }
#         x <- x[order(x[,stcol]),]
#         chrs <- ordered(x[,chrcol], levels=chrorder)
#         # chrs <- ordered(tmp, levels=paste("chr", c(seq(22), "X", "Y"), sep=""))
#         x <- x[order(chrs),]
# }

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
mle_output_file <- "${bam1.baseName}_VS_${bam2.baseName}_mle.bw"

################################################
# DEFAULT PARAMETERS
################################################

remove_anomalies <- FALSE
debug_mode <- TRUE
stebp = 50

################################################
# CORE PROCESSES
################################################

if (debug_mode){
    print(c("ip_file=",ip_file))
    print(c("input_file=",input_file))
    print(c("chromsizes_file=",chromsizes_file))
    print(c("mle_output_file",mle_output_file))
}

## PROCESSING

### Make the list of chromosomes

print("first part")
chromsizes <- fread(chromsizes_file, data.table = FALSE)
#chromsizes <- sortbychr(chromsizes, chrcol="V1", stcol="V2", chrorder=paste("chr", c(seq(19), "X", "Y"), sep=""))
#chrs <- as.character(sub('chr', '', chromsizes[, 1]))
chrs <- chromsizes[, 1]
rownames(chromsizes) <- chrs

### Import the data

ip    <- read.bam.tags(ip_file)
input <- read.bam.tags(input_file)

### Select the informative tags
ip_tags <- ip[['tags']]
input_tags <- input[['tags']]

print("second part")
### Remove the tags with anomalies
if(remove_anomalies==TRUE){
    ip_rm <- remove.local.tag.anomalies(ip_tags[chrs])
    input_rm <- remove.local.tag.anomalies(input_tags[chrs])
}else{
    ip_rm <- ip_tags[chrs]
    input_rm <- input_tags[chrs]
}

print("third part")
### Make the smoothing with the loglikelihood function
print("line 90")
#print(paste0("ip_rm=",head(ip_rm)))
#print(paste0("input_rm=",head(ip_rm)))
# chip_smoothed_mle <- get.smoothed.enrichment.mle(ip_rm, input_rm, tag.shift = 0, background.density.scaling = FALSE)
ip_chrs <- names(ip_rm)[sapply(ip_rm, length) > 10]
input_chrs <- names(input_rm)[sapply(input_rm, length) > 10]
chr2keep <- intersect(ip_chrs, input_chrs)
# length(chr2keep)

chip_smoothed_mle <- get.smoothed.enrichment.mle(ip_rm[chr2keep], input_rm[chr2keep], tag.shift = 0, background.density.scaling = FALSE, step = stebp)

print("line 91")
#chip_smoothed_mle <- get.smoothed.enrichment.mle(ip_tags, input_tags, tag.shift = 0, background.density.scaling = FALSE)

ccl <- sapply(chip_smoothed_mle, nrow)
cpos <- unlist(sapply(chip_smoothed_mle, function(csm) csm\$x))
cend <- cpos + stebp - 1
ccov <- unlist(sapply(chip_smoothed_mle, function(csm) csm\$y))
print("line 97")
gr.mle <- GRanges(seqnames = rep(names(chip_smoothed_mle), ccl), ranges = IRanges(start = cpos, end = cend), strand = Rle('*', sum(ccl)), score = ccov)
seqlengths(gr.mle) <- chromsizes[chr2keep, 2]

#rename score col
col_names <- names(mcols(gr.mle))
#col_names <- sub("\\..*\$", "", col_names)
col_names <- "score"
names(mcols(gr.mle)) <- col_names

print("fouth part")
#export out file
export.bw(gr.mle, mle_output_file)


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

# sink(paste(output_prefix, "R_sessionInfo.log", sep = '.'))
# print(sessionInfo())
# sink()

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
