library(optparse)
options(warn = -1)

# ========================================================================================
# parse command line variables
# ========================================================================================
option_list = list(
    make_option(
        c("-p", "--pop"),
        type    = "character",
        default = NULL, 
        help    = "Name of the population used for building prediction models",
        metavar = "character"
    ),
    make_option(
        c("-a", "--altpop"),
        type    = "character",
        default = NULL, 
        help    = "Name of the population used for validating prediction models",
        metavar = "character"
    ),
    make_option(
        c("-rp", "--rna-pop"),
        type    = "character",
        default = NULL, 
        help    = "RNA data from the population used for building prediction models",
        metavar = "character"
    ),
    make_option(
        c("-ra", "--rna-altpop"),
        type    = "character",
        default = NULL, 
        help    = "RNA data from the population used for validating prediction models",
        metavar = "character"
    ),
    make_option(
        c("-pp", "--predictions-pop"),
        type    = "character",
        default = NULL, 
        help    = "Predictions from the population used for building prediction models",
        metavar = "character"
    ),
    make_option(
        c("-pa", "--predictions-altpop"),
        type    = "character",
        default = NULL, 
        help    = "Predictions from the population used for validating prediction models",
        metavar = "character"
    ),
    make_option(
        c("-o", "--output-prefix"),
        type    = "character",
        default = NULL, 
        help    = "Prefix for output files",
        metavar = "character"
    ),
    make_option(
        c("-k", "--sample-pop-key"),
        type    = "character",
        default = NULL, 
        help    = "Path to two-column tab-separated file with SAMPLE and POP code",
        metavar = "character"
    )
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = TRUE)

pop    = opt$pop
altpop = opt$pop

pop.rna.path  = opt$rna_pop
pop.pred.path = opt$predictions_pop

altpop.rna.path  = opt$rna_altpop
altpop.pred.path = opt$predictions_altpop

output.pfx = opt$output_prefix

key.file.path = opt$sample_pop_key

# =======================================================================================
# load libraries
# =======================================================================================
library(methods)
library(data.table)
library(dplyr)
library(purrr)
library(broom)

# =======================================================================================
# subroutines
# =======================================================================================

# subroutine to compute coefficient of determination from linear model
lmtest = function(x, lm.formula) {
    x.nona = na.omit(x)
    n = dim(x.nona)[1]
    if ( n < 1 ) {
        return(data.table(Gene = x$Gene, R2 = NA, N = 0))
    }
    my.lm = lm(formula(lm.formula), data = x.nona, na.action = na.omit)
    return(data.table(Gene = x$Gene, R2 = summary(my.lm)$r.squared, N = n))
}

# subroutine for extracting both Spearman rho, p-value from correlation test
cortest = function(x, cor.formula) {
    x.nona = na.omit(x)
    n = dim(x.nona)[1]
    if ( n < 1 ) {
        return(data.table(Gene = x$Gene, Correlation = NA, Corr.p.value = NA))
    }
    my.cortest = cor.test(formula(cor.formula), x.nona, method = "spearman", na.action = na.omit)
    return(data.table(Gene = x.nona$Gene, Correlation = my.cortest$estimate, Corr.p.value = my.cortest$p.value))
}

compute.r2 = function(x, lm.formula, from.pop, to.pop){
    r2s = x %>% 
        group_by(Gene, Pop) %>% 
        do(lmtest(., lm.formula)) %>% 
        as.data.table %>% 
        unique
    colnames(r2s) = c("Pop", "Gene", paste0("R2_", from.pop, "_to_", to.pop), paste0("Num_Pred_", from.pop, "_to_", to.pop))
    return(r2s)
}

compute.corr = function(x, corr.formula, from.pop, to.pop){
    corrs = x %>% 
        group_by(Gene, Pop) %>% 
        do(cortest(., corr.formula)) %>% # use previous subroutine to perform correlation test and extract three columns (gene, correlation, p-value)
        as.data.table %>% # cast as data.table to purge duplicate rows
        unique # need this because inelegant subroutine prints repeated rows, 1 per sample instead of 1 per gene group
    colnames(corrs) = c("Pop", "Gene" , paste0("Corr_", from.pop, "_to_", to.pop), paste0("Corr_pval_", from.pop, "_to_", to.pop))
    return(corrs)
}

# ========================================================================================
# load files 
# ========================================================================================

pop.rna  = fread(pop.rna.path)
pop.pred = fread(pop.pred.path)

altpop.rna  = fread(altpop.rna.path)
altpop.pred = fread(altpop.pred.path)

# key.file structure; no header, 1st col Subject ID, 2nd Pop
key.file = fread(key.file.path, header = FALSE)
colnames(key.file) = c("SubjectID", "Pop") 

# melt rna files since they are arranged as matrices
pop.rna.melt    = melt(pop.rna, id.vars = "Gene", variable.name = "SubjectID")
altpop.rna.melt = melt(altpop.rna, id.vars = "Gene", variable.name = "SubjectID")
colnames(pop.rna.melt)[3]    = "Measured_Expr"
colnames(altpop.rna.melt)[3] = "Measured_Expr"

# melt predictions if necessary
if (dim(pop.pred)[2] != 3){
    pop.pred = melt(pop.pred, id.vars = "Gene", variable.name = "SubjectID")
}
colnames(pop.pred)[3] = "Predicted_Expr"

if (dim(altpop.pred)[2] != 3){
    altpop.pred = melt(altpop.pred, id.vars = "Gene", variable.name = "SubjectID")
}
colnames(altpop.pred)[3] = "Predicted_Expr"

# merge predictions and measurements
pop.rnapred    = merge(pop.rna.melt, pop.pred, by = c("Gene", "SubjectID"))
altpop.rnapred = merge(altpop.rna.melt, altpop.pred, by = c("Gene", "SubjectID"))

# merge rnapred data.tables
# we will keep populations separate using key.file later
rnapred = rbind(pop.rnapred, altpop.rnapred)

# add population code
# result is a data.frame, recast it to data.table
#pop.rnapred    = as.data.table(left_join(pop.rnapred, key.file, by = "SubjectID"))
#altpop.rnapred = as.data.table(left_join(altpop.rnapred, key.file, by = "SubjectID"))

# add population code
# use a right_join to ensure that we only save predictions in crosspop scheme
# this discards the samples not in the 89 selected for each pop
# result is a data.frame, so recast it to data.table
rnapred = as.data.table(right_join(rnapred, key.file, by = "SubjectID"))

# compute R2s
lm.formula       = "Predicted_Expr ~ Measured_Expr"
#pop.r2s.topop    = compute.r2(pop.rnapred, lm.formula, pop, pop)
#pop.r2s.toaltpop = compute.r2(altpop.rnapred, lm.formula, pop, altpop)
r2s.allpop = compute.r2(rnapred, lm.formula, pop, "allpop")

# compute correlations
cor.formula        = "~ Predicted_Expr + Measured_Expr"
#pop.corrs.topop    = compute.corr(pop.rnapred, cor.formula, pop, pop)
#pop.corrs.toaltpop = compute.corr(altpop.rnapred, cor.formula, pop, altpop)
corrs.allpop = compute.corr(rnapred, cor.formula, pop, "allpop")

# merge results
#pop.topop.results    = merge(pop.r2s.topop, pop.corrs.topop, by = c("Gene"))
#pop.toaltpop.results = merge(pop.r2s.toaltpop, pop.corrs.toaltpop, by = c("Gene"))
allpop.results = merge(r2s.allpop, corrs.allpop, by = c("Gene", "Pop"))

# save results
#pop.r2.path = paste0(output.pfx, ".", pop, ".predictinto.", pop, ".R2.txt")
#altpop.r2.path = paste0(output.pfx, ".", pop, ".predictinto.", altpop, ".R2.txt")
#fwrite(x = pop.r2s.topop, file = pop.r2.path, quote = FALSE, na = "NA", sep = "\t")
#fwrite(x = pop.r2s.toaltpop, file = altpop.r2.path, quote = FALSE, na = "NA", sep = "\t")
#
#pop.corr.path = paste0(output.pfx, ".", pop, ".predictinto.", pop, ".corr.txt")
#altpop.corr.path = paste0(output.pfx, ".", pop, ".predictinto.", altpop, ".corr.txt")
#fwrite(x = pop.corrs.topop, file = pop.corr.path, quote = FALSE, na = "NA", sep = "\t")
#fwrite(x = pop.corrs.toaltpop, file = altpop.corr.path, quote = FALSE, na = "NA", sep = "\t")

#pop.results.path = paste0(output.pfx, ".", pop, ".predictinto.", pop, ".results.txt")
#altpop.results.path = paste0(output.pfx, ".", pop, ".predictinto.", altpop, ".results.txt")
#fwrite(x = pop.topop.results, file = pop.results.path, quote = FALSE, na = "NA", sep = "\t")
#fwrite(x = pop.toaltpop.results, file = altpop.results.path, quote = FALSE, na = "NA", sep = "\t")

allpop.results.path = paste0(output.pfx, ".", pop, ".predictinto.allpop.results.txt")
fwrite(x = allpop.results, file = allpop.results.path, quote = FALSE, na = "NA", sep = "\t")

# proceed to compute summaries of these results
# start by renaming columns for convenience
colnames(allpop.results) = c("Gene", "Pop", "R2", "Num_Pred", "Corr", "Corr_pval")

# get the predicted genes in common across all pops
genes.in.common = allpop.results %>%
    group_by(Gene) %>%
    select(Gene) %>%
    mutate(n = n()) %>%
    dplyr::filter(., n == 5) %>% ### TODO: determine number of pops programmatically? currently set to 5
    distinct %>%
    select(Gene) %>%
    as.data.table

# save these genes to file
commongenes.path = paste0(output.pfx, ".", pop, ".commongenes.txt")
fwrite(x = genes.in.common, file = commongenes.path, quote = FALSE, sep = "\t", col.names = FALSE)

# get mean R2
# first computes mean R2 for each gene in each pop
# then reduces across genes to get average R2 for pop
mean.r2 = allpop.results %>%
    select(Gene, Pop, R2) %>%
    dplyr::filter(., Gene %in% genes.in.common$Gene) %>%
    group_by(Gene, Pop) %>%
    summarize(R2_mean = mean(R2, na.rm = TRUE)) %>%
    group_by(Pop) %>%
    summarize(Mean = mean(R2_mean, na.rm = TRUE),
              Std_Err = sd(R2_mean, na.rm = TRUE),
              Num_Genes = n()
    ) %>%
    as.data.table

# rename columns to clarify the training population
colnames(mean.r2)[2] = paste0("R2_mean_from_", pop) 
colnames(mean.r2)[3] = paste0("R2_stderr_from_", pop) 

# save to file
allpop.summary.path = paste0(output.pfx, ".", pop, ".predictinto.allpop.summary.txt")
fwrite(x = mean.r2, file = allpop.summary.path, quote = FALSE, na = "NA", sep = "\t")

# repeat search of genes in common across pops
# this time purge those with nonpositive correlation
genes.in.common = allpop.results %>%
    group_by(Gene) %>%
    subset(Corr > 0) %>% # <-- operative difference
    select(Gene) %>%
    mutate(n = n()) %>%
    dplyr::filter(., n == 5) %>%
    distinct %>%
    select(Gene) %>%
    as.data.table

# save these genes to file
commongenes.path = paste0(output.pfx, ".", pop, ".commongenes.poscor.txt")
fwrite(x = genes.in.common, file = commongenes.path, quote = FALSE, sep = "\t", col.names = FALSE)

# get another mean R2
# this one discards negative correlations
# produces more realistic R2s but discards a lot of genes
mean.r2 = allpop.results %>%
    select(Gene, Pop, R2) %>%
    group_by(Gene) %>%
    dplyr::filter(., Gene %in% genes.in.common$Gene) %>%
    group_by(Gene, Pop) %>%
    summarize(R2_mean = mean(R2, na.rm = TRUE)) %>%
    group_by(Pop) %>%
    summarize(Mean = mean(R2_mean, na.rm = TRUE),
              Std_Err = sd(R2_mean, na.rm = TRUE),
              Num_Genes = n()
    ) %>%
    as.data.table

# rename columns to clarify the training population
colnames(mean.r2)[2] = paste0("R2_mean_from_", pop) 
colnames(mean.r2)[3] = paste0("R2_stderr_from_", pop) 

# save to file
allpop.summary.path = paste0(output.pfx, ".", pop, ".predictinto.allpop.summary.poscorr.txt")
fwrite(x = mean.r2, file = allpop.summary.path, quote = FALSE, na = "NA", sep = "\t")
