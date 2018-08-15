#!/usr/bin/env Rscript --vanilla

library(data.table)

# parse command line arguments
args = commandArgs(trailingOnly = TRUE)

# error check
if (length(args) != 3) {
    stop("Usage: Rscript geuvadis_subsample_eur373.R $NUM_RESAMPLE $EUR373_RNA_PATH $EUR89_OUTPUT_DIR")
}

nresample       = as.numeric(args[1])  ## how many subsamples to build?
eur373.rna.path = args[2]  ## where is EUR373 RNA file?
eur89.dir       = args[3]  ## where will EUR89 subsample output go?

# script variables
resample.size = 89
seed = 2018  ## <-- VERY IMPORTANT FOR REPRODUCIBILITY

# set seed
set.seed(seed)  ## <-- ALSO VERY IMPORTANT FOR REPRODUCIBILITY

# load expression data
eur373.rna = fread(eur373.rna.path)

# build each subsample in the same manner
for (i in 1:nresample) {

    # each subsampled set needs its own file paths
    eur89.rna.path        = file.path(eur89.dir, paste0("geuvadis.eur89_", i, ".RPKM.invnorm.txt"))
    eur89.pheno.path      = file.path(eur89.dir, paste0("geuvadis.eur89_", i, ".RPKM.invnorm.pheno"))
    eur89.subjectids.path = file.path(eur89.dir, paste0("geuvadis.eur89_", i, ".sampleids.txt"))

    # this variable will house the "Gene" identifiers and a smattering of sample column names
    # can use this to (consistently) pull data from eur373.rna
    # will also use this to construct files needed for PLINK
    new.colnames = c("Gene", sample(colnames(eur373.rna)[-1], resample.size))

    # subset the data
    eur89.rna = eur373.rna[, new.colnames, with=FALSE ]

    # write this object to file; this is the TXT variant
    fwrite(x = eur89.rna, file = eur89.rna.path, sep = "\t", quote = FALSE)

    # melt and recast the data into a PLINK PHENO format
    # PHENO format requires FID, IID columns as 1st,2nd columns of data frame
    # here, both columns are SubjectID
    eur89.rna.melt = melt(eur89.rna, variable.name = "SubjectID")
    eur89.pheno    = dcast(eur89.rna.melt, SubjectID ~ Gene)
    eur89.pheno    = cbind(eur89.pheno$SubjectID, eur89.pheno)

    # write PHENO variant to file
    # part of PHENO specification is that file DOES NOT have column names
    fwrite(x = eur89.pheno, file = eur89.pheno.path, sep = "\t", quote = FALSE, col.names = FALSE)

    # finally, write subject IDs in PLINK-readable format
    # this is merely SubjectIDs repeated twice
    # as with PHENO, file DOES NOT have column names
    fwrite(x = data.table("A" = eur89.pheno$SubjectID, "B" = eur89.pheno$SubjectID), file = eur89.subjectids.path, sep = "\t", quote = FALSE, col.names = FALSE)
}
