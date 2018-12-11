#!/usr/bin/env Rscript --vanilla

library(data.table)
library(methods)

# parse command line arguments
args     = commandArgs(TRUE)
in.file  = args[1] # "GEUVADIS.chr22.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.rsq_0.8_maf_0.01_hwe_0.00001_geno_0.05.bim"
id.file  = args[2] # "../snp_ids/Phase1.Geuvadis_dbSnp137_idconvert.txt"
out.file = args[3]

# load PLINK BIM
x = fread(in.file)

# load ID conversion file
y = fread(id.file, header = FALSE)

# remap gEUVADIS SNP IDs to rsIDs
rsids = y$V1[y$V2 %in% x$V2]
x$V2  = rsids

# write new BIM to file
fwrite(x = x, file = out.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
