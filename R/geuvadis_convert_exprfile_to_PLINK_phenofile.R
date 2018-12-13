#!/usr/bin/env Rscript --vanilla
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#
# This script converts a file of GEUVADIS inverse rank normalized RPKMs
# and outputs a formatted phenotype file for use with PLINK and GCTA.
#
# Arguments:
# 
# -- $exprfile_in is the input file name of GEUVADIS inverse rank normalized RPKMs
# -- $phenofile_out is the output file name of  
#
# usage:
# > Rscript geuvadis_convert_exprfile_to_PLINK_phenofile.R $exprfile_in $phenofile_out
#
# coded by Kevin L. Keys (2018)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# load library for fread()
library(data.table)

# parse arguments
args = commandArgs(trailingOnly = TRUE)
exprfile.in   = args[1]
phenofile.out = args[2]

# load input
exprfile = fread(exprfile.in)

# transpose data.table
# in process, remove column of gene names
texprfile = data.table(t(exprfile[,-1]))

# get the subject IDs from original data.table
subjectids = colnames(exprfile)[-1]

# append them to left
# must add two columns of subject IDs (FID, IID)
texprfile = data.table(cbind(subjectids, subjectids, texprfile))

# write to file
# note that file does not contain header
fwrite(
    file      = phenofile.out,
    x         = texprfile,
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE,
    sep       = "\t"
)
