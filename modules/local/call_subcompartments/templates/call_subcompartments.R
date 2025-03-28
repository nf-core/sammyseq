#!/usr/bin/env Rscript

# Load libraries
suppressMessages({
  require( data.table )
  require( doParallel )
  require( GenomicRanges )
  require( factoextra )
  require( parallel )
  require( rtracklayer )
  require( ggplot2 )
  require( R.utils )
  require( R.oo )
  require( R.methodsS3 )
  require( maptools )
  require( patchwork )
  require( sp )
  require( dendextend )
  require( Gviz )
  require( CALDER )
  require( ape )
  require( fitdistrplus )
  require( igraph )
  require( Matrix )
  require( rARPACK )
  require( fields )
  require( strawr )
})

args <- commandArgs(trailingOnly = TRUE)

tsv_file <- args[1]

comp_df <- fread(tsv_file, data.table = FALSE)

print(colnames(comp_df))


print(head(comp_df, 3))

output_file <- "test_compartments_output.txt"


