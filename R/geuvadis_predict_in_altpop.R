##!/usr/bin/env Rscript --vanilla

suppressMessages(library(data.table))
suppressMessages(library(dplyr))

# parse command line arguments
args = commandArgs(trailingOnly = TRUE)
weightsfile      = args[1]
geno.altpop.file = args[2]
altpop.pred.file = args[3]
gene             = args[4]

# read prediction weights from file
weights = fread(weightsfile)

# purge missing values
weights = na.omit(weights)

# average across all weights for each SNP
betas = data.table(weights %>% group_by(SNP) %>% summarize(Beta = mean(Beta)))

# load PLINK RAW file from alt pop
geno.altpop = fread(geno.altpop.file)

# PLINK RAW normally adds minor/counting allele to marker names
# we need to purge it
colnames(geno.altpop) = gsub("_.*", "", colnames(geno.altpop))

# subset the genotype and weights data.tables to include just the relevant markers
geno.altpop.sub = subset(geno.altpop, select = c("IID", intersect(colnames(geno.altpop), betas$SNP)))
betas.sub = betas[betas$SNP %in% colnames(geno.altpop.sub),]

# set column order of geno.altpop to match the *row* order of the SNP column in betas.sub 
setorder(betas.sub, SNP)
setcolorder(geno.altpop.sub, c("IID", betas.sub$SNP))

# predictions are now a simple matrix-vector operation
geno.altpop.sub$Prediction_into_Altpop = as.matrix(geno.altpop.sub[,-1]) %*% as.matrix(betas.sub$Beta)

# subset predictions and save to file
altpop.pred = geno.altpop.sub[,.(IID, Prediction_into_Altpop)]
altpop.pred$Gene = gene
altpop.pred = altpop.pred[,.(Gene, IID, Prediction_into_Altpop)]
#colnames(altpop.pred) = c("SubjectID", "Gene", "Prediction_into_Altpop")
fwrite(x = altpop.pred, file = altpop.pred.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t") 
