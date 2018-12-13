#!/usr/bin/env Rscript --vanilla
# ==============================================================================
# coded by Kevin L. Keys (2018)
#
# This script will compile transcriptome imputation results,
# using paired genotype-expression data from the GEUVADIS study.
# This analysis matches sample sizes and training power between EUR and AFR training sets.
#
# The script assumes several completed analyses:
# -- imputation from AFR (YRI89) into all EUR (EUR373)
# -- subsampling of EUR into various (and potentially overlapping) groups of 89 subjects
#
# The subsampling analysis originally performed 100 subsamples.
# Control the number of subsamples analyzed here with the script variable "nsamples".
#
# This script saves results along the way, but will not save the final result.
# Instead, it is printed to the console as two tables.
# The first contains average cross-population imputation R2.
# The second contains standard errors from those mean R2. 
# 
# Usage:
#     Rscript geuvadis_summarize_eur89_subsampling_analysis.R
# ==============================================================================

# ==============================================================================
# load libraries
# ==============================================================================
library(data.table)
library(dplyr)
library(broom)

# ==============================================================================
# function definitions
# ==============================================================================

# subroutine for extracting both Spearman rho, p-value from correlation test
cortest = function(x) {
    my.cortest =  cor.test(~ Predicted_Expr + Measured_Expr, data = x, method = "spearman")
    return(data.table(Gene = x$Gene, Correlation = my.cortest$estimate, Corr.p.value = my.cortest$p.value))
}

# subroutine to wrangle Spearman correlations for each gene in a measurement/prediction sample
# relies on previous cortest() subroutine
compute.correlations = function(df, new.colnames) {
    new.df = df %>%
        group_by(Gene) %>%
        na.omit %>%
        do(cortest(.)) %>% # use previous subroutine to perform correlation test and extract three columns (gene, correlation, p-value)
        as.data.table %>% # cast as data.table to purge duplicate rows
        unique # need this because inelegant subroutine prints repeated rows, 1 per sample instead of 1 per gene group
    colnames(new.df) = new.colnames
    return(new.df)
}

# subroutine to compute R2 for each gene in a measurement/prediction sample
# similar in nature to previous compute.correlations() subroutine
compute.r2 = function(df, new.colnames) {
    new.df = df %>%
        group_by(Gene) %>%
        na.omit %>%
        do(tidy(summary(lm(Predicted_Expr ~ Measured_Expr, data = .))$r.squared)) %>%
        select(Gene, x) %>%
        as.data.table
    colnames(new.df) = new.colnames
    return(new.df)
} 

# ==============================================================================
# script variables
# ==============================================================================

# script constants
nsamples = 10 # how many subsamples did we do?


# ==============================================================================
# script paths and directories 
# ==============================================================================
HOME             = "/netapp/home/kkeys/"
geuvadis.dir     = file.path(HOME, "gala_sage", "rnaseq", "geuvadis")
afr.resultsdir   = file.path(HOME, "gala_sage", "rnaseq", "glmnet", "elasticnet", "yri89")
eur89.resultsdir = file.path(HOME, "gala_sage", "rnaseq", "glmnet", "elasticnet", "eur89")
yri89.rna.path   = file.path(geuvadis.dir, "geuvadis.yri89.RPKM.invnorm.txt")
yri89.pred.path  = file.path(afr.resultsdir, "geuvadis_elasticnet_yri89_predictions.txt")
yri89.into.eur373.path        = file.path(afr.resultsdir, "geuvadis_elasticnet_yri89_predictinto_eur373.txt")
yri89.into.yri89.results.path = file.path(afr.resultsdir, "geuvadis_elasticnet_yri89_crosspop_yri89_prediction_results.txt")

# ==============================================================================
# start script execution
# ==============================================================================

# load data in common between all subsamples
yri89.into.eur373 = fread(yri89.into.eur373.path)
yri89.rna         = fread(yri89.rna.path)
yri89.pred        = fread(yri89.pred.path)

# want easy, memorable column names for AFR -> EUR predictions
colnames(yri89.into.eur373) = c("Gene", "SubjectID", "Predicted_Expr")

# apply usual scheme to RNA measurements and glmnet predictions
# molted data.table of measurements is easier to merge to predictions
yri89.rna.melt   = melt(yri89.rna, id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr")
yri89.into.yri89 = melt(yri89.pred, id.vars = "Gene", variable.name = "SubjectID", value.name = "Predicted_Expr")

# compile prediction results AFR -> AFR
# compute R2 and correlations separately
# merge them and save to file
cat(paste0("\n\tanalyzing predictions for AFR into AFR...\n")) 
yri89.into.yri89 = merge(yri89.into.yri89, yri89.rna.melt, by = c("Gene", "SubjectID"))

yri89.into.yri89.r2.colnames = c("Gene", "R2_yri89_into_yri89")
yri89.into.yri89.r2 = compute.r2(yri89.into.yri89, yri89.into.yri89.r2.colnames)

yri89.into.yri89.corr.colnames = c("Gene", "Correlation_yri89_into_yri89", "Corr_p.value_yri89_into_yri89")
yri89.into.yri89.corr = compute.correlations(yri89.into.yri89, yri89.into.yri89.corr.colnames)

yri89.into.yri89.results = merge(yri89.into.yri89.r2, yri89.into.yri89.corr, by = "Gene")

fwrite(yri89.into.yri89.results, file = yri89.into.yri89.results.path, col.names = TRUE, quote = FALSE, na = "NA")

afr2afr.summary = unname(apply(na.omit(yri89.into.yri89.results)[,-1], 2, function(x) mean(as.numeric(x))))

# clear some memory before running 
yri89.into.yri89      = FALSE
yri89.into.yri89.r2   = FALSE
yri89.into.yri89.corr = FALSE

# create population labels, one for each analyzed subsample
pops = paste0("eur89_", 1:nsamples)

# preallocate data.table and fill results by row 
# columns are: afr2eur (r2, rho, rho p-val), eur2afr (same 3 columns), eur2eur (same cols), and afr2afr (same cols) 
# will eventually fill with results from subsample eur89_1, eur89_2, etc. 
results.per.pop = data.table(
    "Subsample"            = rep("", nsamples),
    "R2_AFR2EUR"           = rep(-Inf, nsamples),
    "Corr_AFR2EUR"         = rep(-Inf, nsamples),
    "Corr_p.value_AFR2EUR" = rep(-Inf, nsamples),
    "R2_EUR2AFR"           = rep(-Inf, nsamples),
    "Corr_EUR2AFR"         = rep(-Inf, nsamples),
    "Corr_p.value_EUR2AFR" = rep(-Inf, nsamples),
    "R2_EUR2EUR"           = rep(-Inf, nsamples),
    "Corr_EUR2EUR"         = rep(-Inf, nsamples),
    "Corr_p.value_EUR2EUR" = rep(-Inf, nsamples),
    "R2_AFR2AFR"           = rep(afr2afr.summary[1], nsamples),
    "Corr_AFR2AFR"         = rep(afr2afr.summary[2], nsamples),
    "Corr_p.value_AFR2AFR" = rep(afr2afr.summary[3], nsamples)
)

### TODO: CREATE FUNCTION FOR PARTS BELOW?
# loop over all subsamples EUR89_1, EUR89_2, ..., EUR89_NSAMPLE
for (i in 1:nsamples) {

    # which pop are we analyzing?
    pop = pops[i]
    cat(paste0("\nCompiling results for subsample ", pop, "...\n"))

    # files paths particular to current subsample 
    pop.dir = file.path(eur89.resultsdir, pop)
    eur89.ids.path = file.path(geuvadis.dir, "eur89", paste0("geuvadis.", pop, ".sampleids.txt"))
    eur89.into.yri89.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_predictinto_yri89.txt"))
    #eur89.into.eur89.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_predictinto_", pop, ".txt"))
    eur89.into.eur89.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_predictions.txt"))
    eur89.rna.path = file.path(geuvadis.dir, "eur89", paste0("geuvadis.", pop, ".RPKM.invnorm.txt"))
    eur89.into.yri89.results.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_crosspop_yri89_prediction_results.txt"))
    eur89.into.eur89.results.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_crosspop_", pop, "_prediction_results.txt"))

    # load data particular to current subsample
    eur89.ids = fread(eur89.ids.path, header = FALSE)[[1]]  ## subj ids are repeated in 2 cols, only need 1 col
    eur89.rna = fread(eur89.rna.path) 
    eur89.into.yri89 = fread(eur89.into.yri89.path)
    eur89.into.eur89.df = fread(eur89.into.eur89.path)

    # melt the predictions
    eur89.into.eur89 = melt(eur89.into.eur89.df, id.vars = "Gene", variable.name = "SubjectID", value.name = "Predicted_Expr")

    # ensure that we know what each column is
    # want to eventually preserve subsample results with suffix EUR89_N
    colnames(eur89.into.eur89) = c("Gene", "SubjectID", "Predicted_Expr")
    colnames(eur89.into.yri89) = c("Gene", "SubjectID", "Predicted_Expr")

    # apply usual scheme to RNA measurements
    # molted data.table of measurements is easier to merge to predictions
    eur89.rna.melt = melt(eur89.rna, id.vars = "Gene", variable.name = "SubjectID", value.name = "Measured_Expr")

    # pull out AFR prediction into this subsample
    yri89.into.eur89 = subset(yri89.into.eur373, subset = SubjectID %in% eur89.ids)

    # merge predictions with measurements
    yri89.into.eur89 = merge(yri89.into.eur89, eur89.rna.melt, by = c("Gene", "SubjectID"), all = TRUE)
    eur89.into.yri89 = merge(eur89.into.yri89, yri89.rna.melt, by = c("Gene", "SubjectID"), all = TRUE)
    eur89.into.eur89 = merge(eur89.into.eur89, eur89.rna.melt, by = c("Gene", "SubjectID"), all = TRUE)

    ###
    ### compute summaries of AFR -> EUR
    ###
    cat(paste0("\n\tanalyzing predictions for AFR into ", pop, "...\n")) 

    # compute the R2
    yri89.into.eur89.r2.colnames = c("Gene", paste0("R2_yri89_into_", pop))
    yri89.into.eur89.r2 = compute.r2(yri89.into.eur89, yri89.into.eur89.r2.colnames)

    # compute the correlations
    yri89.into.eur89.corr.colnames = c("Gene", paste0("Correlation_yri89_into_", pop), paste0("Corr_p.value_yri89_into_", pop))
    yri89.into.eur89.corr = compute.correlations(yri89.into.eur89, yri89.into.eur89.corr.colnames) 

    # merge results
    yri89.into.eur89.results = merge(yri89.into.eur89.r2, yri89.into.eur89.corr, by = "Gene")


    ###
    ### repeat for EUR -> EUR 
    ###
    cat(paste0("\n\tanalyzing predictions for ", pop, " into ", pop, "...\n")) 
    eur89.into.eur89.r2.colnames = c("Gene", paste0("R2_", pop, "_into_", pop))
    eur89.into.eur89.r2 = compute.r2(eur89.into.eur89, eur89.into.eur89.r2.colnames)

    eur89.into.eur89.corr.colnames = c("Gene", paste0("Correlation_", pop, "_into_", pop), paste0("Corr_p.value_", pop, "_into_", pop))
    eur89.into.eur89.corr = compute.correlations(eur89.into.eur89, eur89.into.eur89.corr.colnames) 

    eur89.into.eur89.results = merge(eur89.into.eur89.r2, eur89.into.eur89.corr, by = "Gene")

    fwrite(eur89.into.eur89.results, file = eur89.into.eur89.results.path, col.names = TRUE, quote = FALSE, na = "NA")


    ###
    ### repeat for EUR -> AFR
    ###
    cat(paste0("\n\tanalyzing predictions for ", pop, " into AFR...\n")) 
    eur89.into.yri89.r2.colnames = c("Gene", paste0("R2_", pop, "_into_yri89"))
    eur89.into.yri89.r2 = compute.r2(eur89.into.yri89, eur89.into.yri89.r2.colnames) 

    eur89.into.yri89.corr.colnames = c("Gene", paste0("Correlation_", pop, "_into_yri89"), paste0("Corr_p.value_", pop, "_into_yri89"))
    eur89.into.yri89.corr = compute.correlations(eur89.into.yri89, eur89.into.yri89.corr.colnames) 

    # merge results
    eur89.into.yri89.results = merge(eur89.into.yri89.r2, eur89.into.yri89.corr, by = "Gene", all = TRUE)

    # save those results to file
    fwrite(eur89.into.yri89.results, file = eur89.into.yri89.results.path, col.names = TRUE, quote = FALSE, na = "NA")

    # compute a summary for whole subsample 
    # add results to results.per.pop 
    eur2afr.summary = unname(apply(na.omit(eur89.into.yri89.results)[,-1], 2, function(x) mean(as.numeric(x))))
    afr2eur.summary = unname(apply(na.omit(yri89.into.eur89.results)[,-1], 2, function(x) mean(as.numeric(x))))
    eur2eur.summary = unname(apply(na.omit(eur89.into.eur89.results)[,-1], 2, function(x) mean(as.numeric(x))))
    set(results.per.pop, i, 1L,  pop)
    set(results.per.pop, i, 2L,  eur2afr.summary[1])
    set(results.per.pop, i, 3L,  eur2afr.summary[2])
    set(results.per.pop, i, 4L,  eur2afr.summary[3])
    set(results.per.pop, i, 5L,  afr2eur.summary[1])
    set(results.per.pop, i, 6L,  afr2eur.summary[2])
    set(results.per.pop, i, 7L,  afr2eur.summary[3])
    set(results.per.pop, i, 8L,  eur2eur.summary[1])
    set(results.per.pop, i, 9L,  eur2eur.summary[2])
    set(results.per.pop, i, 10L, eur2eur.summary[3])

    # recover some memory before starting new iteration of loop
    eur89.into.yri89.r2   = FALSE
    eur89.into.yri89.corr = FALSE
    eur89.into.eur89.r2   = FALSE
    eur89.into.eur89.corr = FALSE
    yri89.into.eur89.r2   = FALSE
    yri89.into.eur89.corr = FALSE
    yri89.into.eur89      = FALSE 
    eur89.into.yri89      = FALSE
    gc()
    
}

# compute mean of means for each train/test scenario
results.per.pop.summary.mean.r2 = matrix(Inf, 2, 2)
colnames(results.per.pop.summary.mean.r2) = c("AFR", "EUR")
rownames(results.per.pop.summary.mean.r2) = c("AFR", "EUR")
results.per.pop.summary.mean.r2[1,1] = mean(results.per.pop$R2_AFR2AFR, na.rm = TRUE) 
results.per.pop.summary.mean.r2[1,2] = mean(results.per.pop$R2_EUR2AFR, na.rm = TRUE)  
results.per.pop.summary.mean.r2[2,1] = mean(results.per.pop$R2_AFR2EUR, na.rm = TRUE) 
results.per.pop.summary.mean.r2[2,2] = mean(results.per.pop$R2_EUR2EUR, na.rm = TRUE) 

# compute standard errors of those means
results.per.pop.summary.sd.r2 = matrix(Inf, 2, 2)
colnames(results.per.pop.summary.sd.r2) = c("AFR", "EUR")
rownames(results.per.pop.summary.sd.r2) = c("AFR", "EUR")
results.per.pop.summary.sd.r2[1,1] = sd(results.per.pop$R2_AFR2AFR, na.rm = TRUE) 
results.per.pop.summary.sd.r2[1,2] = sd(results.per.pop$R2_EUR2AFR, na.rm = TRUE)  
results.per.pop.summary.sd.r2[2,1] = sd(results.per.pop$R2_AFR2EUR, na.rm = TRUE) 
results.per.pop.summary.sd.r2[2,2] = sd(results.per.pop$R2_EUR2EUR, na.rm = TRUE) 

# print final tables
cat("\n\nMean R2\n")
print(results.per.pop.summary.mean.r2)

cat("\n\nStdErr R2\n")
print(results.per.pop.summary.sd.r2)

# compute mean of means for each train/test scenario
results.per.pop.summary.mean.corr = matrix(Inf, 2, 2)
colnames(results.per.pop.summary.mean.corr) = c("AFR", "EUR")
rownames(results.per.pop.summary.mean.corr) = c("AFR", "EUR")
results.per.pop.summary.mean.corr[1,1] = mean(results.per.pop$Corr_AFR2AFR, na.rm = TRUE) 
results.per.pop.summary.mean.corr[1,2] = mean(results.per.pop$Corr_EUR2AFR, na.rm = TRUE)  
results.per.pop.summary.mean.corr[2,1] = mean(results.per.pop$Corr_AFR2EUR, na.rm = TRUE) 
results.per.pop.summary.mean.corr[2,2] = mean(results.per.pop$Corr_EUR2EUR, na.rm = TRUE) 

# compute standard errors of those means
results.per.pop.summary.sd.corr = matrix(Inf, 2, 2)
colnames(results.per.pop.summary.sd.corr) = c("AFR", "EUR")
rownames(results.per.pop.summary.sd.corr) = c("AFR", "EUR")
results.per.pop.summary.sd.corr[1,1] = sd(results.per.pop$Corr_AFR2AFR, na.rm = TRUE) 
results.per.pop.summary.sd.corr[1,2] = sd(results.per.pop$Corr_EUR2AFR, na.rm = TRUE)  
results.per.pop.summary.sd.corr[2,1] = sd(results.per.pop$Corr_AFR2EUR, na.rm = TRUE) 
results.per.pop.summary.sd.corr[2,2] = sd(results.per.pop$Corr_EUR2EUR, na.rm = TRUE) 

# print final tables
cat("\n\nMean Spearman Correlation\n")
print(results.per.pop.summary.mean.corr)

cat("\n\nStdErr Spearman Correlation\n")
print(results.per.pop.summary.sd.corr)
