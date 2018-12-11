#!/usr/bin/env Rscript --vanilla
#
# this script subsets RNA-Seq data from gEUVADIS for use with training in glmnet
# the data are pre-downloaded from https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/
# the output of this script consists of six files, three for Europeans (EUR) and two for Yorubans (YRI)
# for each ethnicity, the script produces:
#     (1) a matrix of inverse-normal transformed rank-normalized RPKMs
#     (2) a transposed copy of the matrix from (1)
#     (3) a list of sample IDs for matrices (1) and (2)
#
# This script does not transform the genes to normality. Per the GEUVADIS README
#
#     November 5, 2013 update: The file GD462.GeneQuantRPKM.50FN.samplename.resk10.norm.txt.gz that had the normalization as above
#     PLUS an additional transformation of each gene's values to standard normal has been replaced by GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
#
# https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-3/GeuvadisRNASeqAnalysisFiles_README.txt
#
# coded by Kevin L. Keys (2018)

suppressMessages(library(data.table))
suppressMessages(library(methods))
suppressMessages(library(optparse))

## input file paths
#all.rnaseq.data.file = "~/gala_sage/rnaseq/geuvadis/rnaseq/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt"
#sample.ids.file      = "~/gala_sage/rnaseq/geuvadis/genotypes/sample_ids/gevuadis.465.sample.ids.txt"
#ryan.file            = "expressions.genes.yri.matrix.eqtl.txt"
#
## output file paths
#eur.out.file  = "geuvadis.eur373.RPKM.invnorm.txt"
#yri.out.file  = "geuvadis.yri89.RPKM.invnorm.txt"
#teur.out.file = "geuvadis.eur373.RPKM.invnom.transposed.txt"
#teur.out.file = "geuvadis.eur373.RPKM.invnom.transposed.txt"
#eur.ids.file  = "geuvadis.eur373.ids.txt"
#yri.ids.file  = "geuvadis.yri373.ids.txt"

# parse command line variables
option_list = list(
    make_option(
        c("-r", "--rnaseq-data-file"),
        type    = "character",
        default = NULL,
        help    = "The data frame containing all results to analyze.",
        metavar = "character"
    ),
    make_option(
        c("-s", "--sample-ids-file"),
        type    = "character",
        default = NULL,
        help    = "PLINK-format sample IDs file (two columns, one ID per line, repeated in each column).",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-directory"),
        type    = "character",
        default = NULL,
        help    = "Directory for saving output files",
        metavar = "character"
    ),
    make_option(
        c("-eo", "--EUR-out-file"),
        type    = "character",
        default = "geuvadis.eur373.RPKM.invnorm.txt",
        help    = "Output file for parsed EUR RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-yo", "--YRI-out-file"),
        type    = "character",
        default = "geuvadis.yri89.RPKM.invnorm.txt",
        help    = "Output file for parsed YRI RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-teo", "--transposed-EUR-out-file"),
        type    = "character",
        default = "geuvadis.eur373.RPKM.invnorm.transposed.txt",
        help    = "Transposed output file for parsed EUR RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-yo", "--transposed-YRI-out-file"),
        type    = "character",
        default = "geuvadis.yri89.RPKM.invnorm.transposed.txt",
        help    = "Transposed output file for parsed YRI RNA-Seq data [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-ei", "--EUR-IDs-file"),
        type    = "character",
        default = "geuvadis.eur373.ids.txt",
        help    = "Output file for parsed EUR sample IDs [default = %default].",
        metavar = "character"
    ),
    make_option(
        c("-yi", "--YRI-IDs-file"),
        type    = "character",
        default = "geuvadis.yri89.ids.txt",
        help    = "Output file for parsed EUR sample IDs [default = %default].",
        metavar = "character"
    )
)


 opt_parser = OptionParser(option_list = option_list)
 opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

 cat("Parsed options:\n\n")
 print(opt)

# input file paths
all.rnaseq.data.file = opt$rnaseq_data_file 
sample.ids.file      = opt$sample_ids_file 
ryan.file            = args[3]
eur.out.file  = file.path(opt$output_directory, opt$EUR_out_file) 
yri.out.file  = file.path(opt$output_directory, opt$YRI_out_file) 
teur.out.file = file.path(opt$output_directory, opt$transposed_EUR_out_file) 
tyri.out.file = file.path(opt$output_directory, opt$transposed_YRI_out_file) 
eur.ids.file  = file.path(opt$output_directory, opt$EUR_IDs_file) 
yri.ids.file  = file.path(opt$output_directory, opt$YRI_IDs_file) 


# load from file
rnaseq.data = fread(all.rnaseq.data.file)     # load normalized RPKMs for gEUVADIS
ids  = fread(sample.ids.file, header = FALSE)   # load sample IDs
ryan = fread(ryan.file)                         # Ryan Hernandez and Jimmie Ye have done YRI normalization before, we want just the IDs

# discard duplicate gene ID and chr, position info
# rename gene ID to just "Gene"
# also trim any transcript number
rnaseq.data = rnaseq.data[,-c(2:4)]
colnames(rnaseq.data)[1] = "Gene"
rnaseq.data$Gene = strtrim(rnaseq.data$Gene, 15)

# rename gene ID to just "Gene"
colnames(ryan)[1] = "Gene"

# subset the RNA-Seq data
yri = rnaseq.data[,colnames(rnaseq.data) %in% colnames(ryan), with = FALSE]
eur = rnaseq.data[,!(colnames(rnaseq.data) %in% colnames(ryan)), with = FALSE]

# put "Gene" column in eur
eur$Gene = rnaseq.data$Gene

# order eur and yri by gene ID
setorder(x=yri, Gene, na.last = T)
setorder(x=eur, Gene, na.last = T)

# reorder columns of eur, putting "Gene" first
eur = eur[,c(ncol(eur),1:(ncol(eur)-1)), with = F]

# write subsetted data.tables to file
fwrite(x=eur, file = eur.out.file, col.names=T, quote=F)
fwrite(x=yri, file = yri.out.file, col.names=T, quote=F)

# also write transposed data files
# involves transposing data tables and adding back IDs
# do first for eur and then for yri
teur = data.table(t(eur[,-1]))
colnames(teur) = eur$Gene
teur$IID = colnames(eur)[-1]
setkey(teur, IID)
setorder(teur, IID)
teur = teur[,c(ncol(teur),1:(ncol(teur)-1)), with=FALSE]
fwrite(x = teur, file = teur.out.file, quote = FALSE)

tyri = data.table(t(yri[,-1]))
colnames(tyri) = yri$Gene
tyri$IID = colnames(yri)[-1]
setkey(tyri, IID)
setorder(tyri, IID)
tyri = tyri[,c(ncol(tyri),1:(ncol(tyri)-1)), with=FALSE]
fwrite(x = tyri, file = tyri.out.file, quote = FALSE)

# also save separate lists of EUR, YRI sample IDs
# this will come in handy when subsetting genotypes
eur.ids = colnames(eur)[-1]
yri.ids = colnames(yri)[-1]
fwrite(x = eur.ids, file = eur.ids.file, col.names = FALSE, quote = FALSE)
fwrite(x = yri.ids, file = yri.ids.file, col.names = FALSE, quote = FALSE)
