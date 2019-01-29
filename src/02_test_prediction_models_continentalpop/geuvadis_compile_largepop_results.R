#!/usr/bin/env Rscript --vanilla

# ==========================================================================================
# load libraries
# ==========================================================================================
suppressMessages(library(methods))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(broom))
suppressMessages(library(optparse))

# ==========================================================================================
# subroutines
# ==========================================================================================
compute.r2.corr.onegene = function(x){
    the.fit = summary(lm(Predicted_Expr ~ Measured_Expr, data = x)) 
    my.cortest =  cor.test(~ Predicted_Expr + Measured_Expr, data = x, method = "spearman")
    my.output = data.frame(
        "Correlation" = my.cortest$estimate,
        "Corr_pval" = my.cortest$p.value,
        "R2" = the.fit$r.squared,
        "R2_pval" = the.fit$coefficients[2,4]
    )
    return(my.output)
}

compute.r2.corr = function(df){
    new.df = df %>% 
        na.omit %>% 
        dplyr::group_by(Train_Pop, Test_Pop, Gene) %>% 
        do(compute.r2.corr.onegene(.)) %>% 
        as.data.table
    return(new.df)
}

# ==========================================================================================
# parse command line arguments 
# ==========================================================================================

option_list = list(
    make_option(
        c("-a", "--EUR373-to-EUR373"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR373 to all Europeans",
        metavar = "character"
    ),
    make_option(
        c("-b", "--EUR373-to-AFR"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR373 to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-c", "--EUR278-to-EUR278"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to non-Finnish Europeans",
        metavar = "character"
    ),
    make_option(
        c("-d", "--EUR278-to-FIN"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to Finns",
        metavar = "character"
    ),
    make_option(
        c("-e", "--EUR278-to-AFR"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-f", "--AFR-to-EUR373"),
        type    = "character",
        default = NULL,
        help    = "Predictions from EUR278 to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-g", "--AFR-to-AFR"),
        type    = "character",
        default = NULL,
        help    = "Predictions from AFR to the Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-i", "--output-predictions"),
        type    = "character",
        default = NULL,
        help    = "Path to compiled data table of crosspopulation predictions",
        metavar = "character"
    ),
    make_option(
        c("-i", "--output-results"),
        type    = "character",
        default = NULL,
        help    = "Path to compiled data table of prediction results",
        metavar = "character"
    ),
    make_option(
        c("-j", "--output-r2"),
        type    = "character",
        default = NULL,
        help    = "Path to compiled R2 results",
        metavar = "character"
    ),
    make_option(
        c("-k", "--EUR278-ids"),
        type    = "character",
        default = NULL,
        help    = "PLINK-formatted sample IDs for EUR278",
        metavar = "character"
    ),
    make_option(
        c("-l", "--FIN-ids"),
        type    = "character",
        default = NULL,
        help    = "PLINK-formatted sample IDs for FIN",
        metavar = "character"
    ),
    make_option(
        c("-m", "--EUR-RNA"),
        type    = "character",
        default = NULL,
        help    = "Matrix of gene expression measures for Europeans",
        metavar = "character"
    ),
    make_option(
        c("-n", "--AFR-RNA"),
        type    = "character",
        default = NULL,
        help    = "Matrix of gene expression measures for Yoruba",
        metavar = "character"
    ),
    make_option(
        c("-o", "--common-genes"),
        type    = "character",
        default = NULL,
        help    = "Path to save list of genes in common to all train-test scenarios",
        metavar = "character"
    ),
    make_option(
        c("-p", "--common-genes-poscorr"),
        type    = "character",
        default = NULL,
        help    = "Path to save list of genes with positive correlation between data and predictions in all train-test scenarios",
        metavar = "character"
    )
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

cat("Parsed options:\n")
print(opt)


eur373.to.eur373.path = opt$EUR373_to_EUR373
eur373.to.afr.path    = opt$EUR373_to_AFR
eur278.to.eur278.path = opt$EUR278_to_EUR278
eur278.to.fin.path    = opt$EUR278_to_FIN
eur278.to.afr.path    = opt$EUR278_to_AFR
afr.to.eur373.path    = opt$AFR_to_EUR373
#afr.to.eur278.path    = opt$
#afr.to.fin.path       = opt$
eur278.ids.path       = opt$EUR278_ids
fin.ids.path          = opt$FIN_ids
afr.to.afr.path       = opt$AFR_to_AFR
output.path.predictions = opt$output_predictions
output.path.results   = opt$output_results
output.path.r2        = opt$output_r2
eur.rna.path          = opt$EUR_RNA
afr.rna.path          = opt$AFR_RNA
common.genes.path     = opt$common_genes
commongenes.poscorr.path = opt$common_genes_poscorr

# load data frames
eur373.to.eur373 = fread(eur373.to.eur373.path, header = TRUE)
eur373.to.afr    = fread(eur373.to.afr.path, header = TRUE)
eur278.to.eur278 = fread(eur278.to.eur278.path, header = TRUE)
eur278.to.fin    = fread(eur278.to.fin.path, header = TRUE)
eur278.to.afr    = fread(eur278.to.afr.path, header = TRUE)
afr.to.eur373    = fread(afr.to.eur373.path, header = TRUE)
#afr.to.eur278    = fread(afr.to.eur278.path, header = TRUE)
#afr.to.fin       = fread(afr.to.fin.path, header = TRUE)
afr.to.afr       = fread(afr.to.afr.path, header = TRUE)

eur278.ids = fread(eur278.ids.path, header = FALSE)
fin.ids    = fread(fin.ids.path , header = FALSE) 

afr.rna = fread(afr.rna.path, header = TRUE)
eur.rna = fread(eur.rna.path, header = TRUE)

# subset eur278, fin from AFR predictions
afr.to.eur278 = afr.to.eur373 %>% dplyr::filter(SubjectID %in% eur278.ids[[1]]) %>% as.data.table 
afr.to.fin    = afr.to.eur373 %>% dplyr::filter(SubjectID %in% fin.ids[[1]]) %>% as.data.table 

# make char vectors of train/test orders
#pops = c("EUR373_to_EUR373", "EUR373_to_AFR", "EUR278_to_EUR278", "EUR278_to_FIN", "EUR278_to_AFR", "AFR_to_EUR373", "AFR_to_EUR278", "AFR_to_FIN", "AFR_to_AFR")
train.pops = c("EUR373", "EUR373", "EUR278", "EUR278", "EUR278", "AFR", "AFR", "AFR", "AFR")
test.pops = c("EUR373", "AFR", "EUR278", "FIN", "AFR", "EUR373", "EUR278", "FIN", "AFR")

# make list of data frames in same order
datatables = list(eur373.to.eur373, eur373.to.afr, eur278.to.eur278, eur278.to.fin, eur278.to.afr, afr.to.eur373, afr.to.eur278, afr.to.fin, afr.to.afr)

# add pop labels
n = length(train.pops)
for (i in 1:n) { 
    cat("parsing table with train.pop = ", train.pops[i], " and test.pop = ", test.pops[i], "\n")
    #datatables[[i]]$Train_Test = pops[i]
    if ( ncol(datatables[[i]]) != 3) {
        datatables[[i]] = melt(datatables[[i]], id.vars = "Gene", variable.name = "SubjectID", value.name = "Predicted_Expr")
    } else {
        colnames(datatables[[i]]) = c("Gene", "SubjectID", "Predicted_Expr")
    }
    datatables[[i]]$Train_Pop = train.pops[i]
    datatables[[i]]$Test_Pop  = test.pops[i]
    #datatables[[i]]$Spearman.rho = as.numeric(datatables[[i]]$Spearman.rho)
    #setkey(datatables[[i]], Gene, P.value, R2, Spearman.rho, Train_Test) 
    #setkey(datatables[[i]], Gene, P.value, R2, Spearman.rho, Train_Pop, Test_Pop) 
    setkey(datatables[[i]], Gene, SubjectID, Predicted_Expr, Train_Pop, Test_Pop) 
}

# merge data frames together and save to file
geuvadis.predictions = Reduce(function(x,y) merge(x,y,all=TRUE), datatables)
fwrite(x = geuvadis.predictions, file = output.path.predictions, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t") 

# prepare data frame for measured gene expression
afr.rna = melt(afr.rna, id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr") 
eur.rna = melt(eur.rna, id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr") 
rna = rbindlist(list(afr.rna, eur.rna))

# get genes with results in all populations
transcripts = geuvadis.predictions %>%
    na.omit %>%
    count(Gene, Train_Pop, Test_Pop) %>%
    select(Gene) %>%
    count(Gene) %>%
    dplyr::filter(n == n) %>%
    select(Gene) %>%
    as.data.table

# save common transcript IDs to file
fwrite(x = transcripts, file = common.genes.path, row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

# merge predictions + measurements
geuvadis.rnapred = merge(geuvadis.predictions, rna, by = c("Gene", "SubjectID"))

# first compute correlations
geuvadis.results = compute.r2.corr(geuvadis.rnapred)

# use "transcripts" as a filtering variable
#r2 = geuvadis.results %>%
#    dplyr::filter(
#        Gene %in% transcripts[[1]] & 
#        Spearman.rho > 0
#    ) %>%
#    group_by(Train_Test) %>%
#    summarize(r2 = mean(R2, na.rm = T), corr = mean(Spearman.rho, na.rm = T), ntranscripts = length(transcripts[[1]])) %>%
#    as.data.table
#
#r2 = geuvadis.results %>%
#    dplyr::filter(
#        Gene %in% transcripts[[1]] & 
#        Spearman.rho > 0
#    ) %>%
#    group_by(Train_Pop, Test_Pop) %>%
#    summarize(r2 = mean(R2, na.rm = T), corr = mean(Spearman.rho, na.rm = T), ntranscripts = length(transcripts[[1]])) %>%
#    as.data.table

# get transcripts with positive correlation
transcripts.poscorr = geuvadis.results %>%
    na.omit %>%
    dplyr::filter(Correlation > 0) %>%
    count(Gene, Train_Pop, Test_Pop) %>%
    select(Gene) %>%
    count(Gene) %>%
    dplyr::filter(n == n) %>%
    select(Gene) %>%
    as.data.table

# save list of transcripts to file
fwrite(x = transcripts.poscorr, file = common.genes.poscorr.path, row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

r2 = geuvadis.results %>%
    dplyr::filter(Gene %in% transcripts.poscorr[[1]]) %>%
    group_by(Train_Pop, Test_Pop) %>%
    summarize(r2 = mean(R2, na.rm = TRUE), corr = mean(Spearman.rho, na.rm = TRUE), ntranscripts = length(transcripts[[1]])) %>%
    as.data.table

fwrite(x = r2, file = output.path.r2, row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t") 
