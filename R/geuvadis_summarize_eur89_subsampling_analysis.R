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
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(broom))

# ==============================================================================
# function definitions
# ==============================================================================

# subroutine for extracting both Spearman rho, p-value from correlation test
cortest = function(x) {
    my.cortest =  cor.test(~ Predicted_Expr + Measured_Expr, data = x, method = "spearman")
    return(data.table(Gene = x$Gene, Correlation = my.cortest$estimate, Corr.p.value = my.cortest$p.value))
}

compute.r2.corr.onegene = function(x){
    the.fit = summary(lm(Predicted_Expr ~ Measured_Expr, data = x))
    my.cortest =  cor.test(~ Predicted_Expr + Measured_Expr, data = x, method = "spearman")
    my.output = data.frame("Correlation" = my.cortest$estimate, "Corr_pval" = my.cortest$p.value, "R2" = the.fit$r.squared, "R2_pval" = the.fit$coefficients[2,4])
    return(my.output)
}

compute.r2.corr = function(df){
    new.df = df %>%
        na.omit %>%
        dplyr::group_by(Gene) %>%
        do(compute.r2.corr.onegene(.)) %>%
        as.data.table
    return(new.df)
}

summarize.results = function(df) {
    new.df = df %>%
        dplyr::filter(Correlation > 0) %>% 
        summarize(R2 = mean(R2, na.rm = TRUE), Correlation = mean(Correlation, na.rm = TRUE)) %>%
        as.data.table %>%
        unlist
    return(new.df)
}

save.results = function(df, filename, train.pop, test.pop, subsample.num) {
    new.df       = df
    df$Train_Pop = train.pop
    df$Test_Pop  = test.pop
    df$Subsample = subsample.num
    fwrite(new.df, file = filename, row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA", sep = "\t")
    return()
}

summarize.results.bypop.r2 = function(x, train.pop, test.pop) {
    new.results = x %>%
        dplyr::filter(Train_Pop == train.pop & Test_Pop == test.pop) %>%
        summarize(R2_mean = mean(R2, na.rm = TRUE), R2_SD = sd(R2, na.rm = TRUE), Nmodels = n()) %>%
        unlist
    return(new.results)
}

summarize.results.bypop.corr = function(x, train.pop, test.pop) {
    new.results = x %>%
        dplyr::filter(Train_Pop == train.pop & Test_Pop == test.pop) %>%
        summarize(Corr_mean = mean(Correlation, na.rm = TRUE), Corr_SD = sd(Correlation, na.rm = TRUE), Nmodels = n()) %>%
        unlist
    return(new.results)
}

# ==============================================================================
# script variables
# ==============================================================================

# script constants
nsamples = 10 # how many subsamples did we do?


# ==============================================================================
# script paths and directories 
# ==============================================================================
HOME             = "/netapp/home/kkeys"
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
# compute R2 and correlations in one pass 
cat("\n\tanalyzing predictions for AFR into AFR...\n")
yri89.into.yri89 = merge(yri89.into.yri89, yri89.rna.melt, by = c("Gene", "SubjectID"))

yri89.into.yri89.results = compute.r2.corr(yri89.into.yri89)
afr2afr.summary = summarize.results(yri89.into.yri89.results)
save.results(yri89.into.yri89.results,  yri89.into.yri89.results.path, "AFR", "AFR", 0)

# clear some memory before running 
yri89.into.yri89      = FALSE

# create population labels, one for each analyzed subsample
pops = paste0("eur89_", 1:nsamples)


# preallocate data.table for results
# need 3 rows per subsample: EUR2EUR, EUR2AFR, AFR2EUR
# will add AFR2AFR results last
nresults.rows = 3*nsamples
results.per.pop = data.table(
    "R2"          = rep(-Inf, nresults.rows), 
    "Correlation" = rep(-Inf, nresults.rows), 
    "Train_Pop"   = rep("",   nresults.rows), 
    "Test_Pop"    = rep("",   nresults.rows), 
    "Subsample"   = rep(-1,   nresults.rows) 
)


### TODO: CREATE FUNCTION FOR PARTS BELOW?
# loop over all subsamples EUR89_1, EUR89_2, ..., EUR89_NSAMPLE
for (i in 1:nsamples) {

    # which pop are we analyzing?
    pop = pops[i]
    cat("\nCompiling results for subsample ", pop, "...\n")

    # files paths particular to current subsample 
    pop.dir = file.path(eur89.resultsdir, pop)
    eur89.ids.path = file.path(geuvadis.dir, "eur89", paste0("geuvadis.", pop, ".sampleids.txt"))
    eur89.into.yri89.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_predictinto_yri89.txt"))
    eur89.into.eur89.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_predictions.txt"))
    eur89.rna.path = file.path(geuvadis.dir, "eur89", paste0("geuvadis.", pop, ".RPKM.invnorm.txt"))
    eur89.into.yri89.results.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_crosspop_yri89_prediction_results.txt"))
    eur89.into.eur89.results.path = file.path(pop.dir, paste0("geuvadis_elasticnet_", pop, "_crosspop_", pop, "_prediction_results.txt"))
    yri89.into.eur89.results.path = file.path(pop.dir, paste0("geuvadis_elasticnet_yri89_crosspop_", pop, "_prediction_results.txt"))

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

    # no longer need RNA data, recover memory
    eur89.rna.melt = FALSE
    gc()

    ###
    ### compute summaries of AFR -> EUR
    ###
    cat("\n\tanalyzing predictions for AFR into ", pop, "...\n")
    yri89.into.eur89.results = compute.r2.corr(yri89.into.eur89)
    afr2eur.summary = summarize.results(yri89.into.eur89.results)
    save.results(yri89.into.eur89.results,  yri89.into.eur89.results.path, "AFR", "EUR", i)

    # no longer need predictions, so destroy and recover memory 
    yri89.into.eur89 = FALSE
    gc()

    ###
    ### repeat for EUR -> EUR 
    ###
    cat("\n\tanalyzing predictions for ", pop, " into ", pop, "...\n")
    eur89.into.eur89.results = compute.r2.corr(eur89.into.eur89)
    eur2eur.summary = summarize.results(eur89.into.eur89.results)
    save.results(eur89.into.eur89.results,  eur89.into.eur89.results.path, "EUR", "EUR", i)
    eur89.into.eur89 = FALSE
    gc()


    ###
    ### repeat for EUR -> AFR
    ###
    cat("\n\tanalyzing predictions for ", pop, " into AFR...\n")
    eur89.into.yri89.results = compute.r2.corr(eur89.into.yri89)
    eur2afr.summary = summarize.results(eur89.into.yri89.results)
    save.results(eur89.into.yri89.results,  eur89.into.yri89.results.path, "EUR", "AFR", i)
    eur89.into.yri89 = FALSE
    gc()


    # compute a summary for whole subsample 
    # add results to results.per.pop 
    row1 = as.integer(i)
    row2 = as.integer(2*i)
    row3 = as.integer(3*i)

    set(results.per.pop, row1, 1L,  afr2eur.summary["R2"])
    set(results.per.pop, row1, 2L,  afr2eur.summary["Correlation"])
    set(results.per.pop, row1, 3L,  "AFR")
    set(results.per.pop, row1, 4L,  "EUR")
    set(results.per.pop, row1, 5L,  i)

    set(results.per.pop, row2, 1L,  eur2eur.summary["R2"])
    set(results.per.pop, row2, 2L,  eur2eur.summary["Correlation"])
    set(results.per.pop, row2, 3L,  "EUR")
    set(results.per.pop, row2, 4L,  "EUR")
    set(results.per.pop, row2, 5L,  i)

    set(results.per.pop, row3, 1L,  eur2afr.summary["R2"])
    set(results.per.pop, row3, 2L,  eur2afr.summary["Correlation"])
    set(results.per.pop, row3, 3L,  "EUR")
    set(results.per.pop, row3, 4L,  "AFR")
    set(results.per.pop, row3, 5L,  i)

    # recover any available memory before starting new iteration of loop
    gc()
}

# save intermediate results
fwrite(x = results.per.pop, file = file.path(eur89.resultsdir, "geuvadis.subsampling.results.bypop.txt"), sep = "\t", quote = FALSE, na = "NA")

cat("creating cross-population summaries of results...\n") 
eur2afr.r2 = summarize.results.bypop.r2(results.per.pop, "EUR", "AFR")
eur2eur.r2 = summarize.results.bypop.r2(results.per.pop, "EUR", "EUR")
afr2eur.r2 = summarize.results.bypop.r2(results.per.pop, "AFR", "EUR")

# expected output:
#> r2.results
#  Train_Pop Test_Pop Measure_Type Measure_Value
#1       EUR      AFR      R2_Mean            -1
#2       EUR      AFR        R2_SD            -1
#3       EUR      EUR      R2_Mean            -1
#4       EUR      EUR        R2_SD            -1
#5       AFR      EUR      R2_Mean            -1
#6       AFR      EUR        R2_SD            -1
#7       AFR      AFR      R2_Mean            -1
#8       AFR      AFR        R2_SD            -1
r2.results = data.frame(
    "Train_Pop" = c(rep("EUR", 4), rep("AFR", 4)),
    "Test_Pop"  = c(rep("AFR", 2), rep("EUR", 4), rep("AFR", 2)),
    "Measure_Type" = rep(c("R2_Mean", "R2_SD"), 4),
    "Measure_Value" = rep(-1, 8)
)
r2.results$Measure_Value[1] = eur2afr.r2["R2_mean"] 
r2.results$Measure_Value[2] = eur2afr.r2["R2_SD"]
r2.results$Measure_Value[3] = eur2eur.r2["R2_mean"]
r2.results$Measure_Value[4] = eur2eur.r2["R2_SD"]
r2.results$Measure_Value[5] = afr2eur.r2["R2_mean"]
r2.results$Measure_Value[6] = afr2eur.r2["R2_SD"]
r2.results$Measure_Value[7] = afr2afr.summary["R2"]
r2.results$Measure_Value[8] = 0

# repeat for correlations
# will make another data.frame for those results
# then slap R2, Correlation results together and save to file

eur2afr.corr = summarize.results.bypop.corr(results.per.pop, "EUR", "AFR")
eur2eur.corr = summarize.results.bypop.corr(results.per.pop, "EUR", "EUR")
afr2eur.corr = summarize.results.bypop.corr(results.per.pop, "AFR", "EUR")

# same as r2.results but "Measure_Type" alternates b/w "Corr_Mean" and "Corr_SD"
corr.results = data.frame(
    "Train_Pop" = c(rep("EUR", 4), rep("AFR", 4)),
    "Test_Pop"  = c(rep("AFR", 2), rep("EUR", 4), rep("AFR", 2)),
    "Measure_Type" = rep(c("Corr_Mean", "Corr_SD"), 4),
    "Measure_Value" = rep(-1, 8)
)
corr.results$Measure_Value[1] = eur2afr.corr["Corr_mean"] 
corr.results$Measure_Value[2] = eur2afr.corr["Corr_SD"]
corr.results$Measure_Value[3] = eur2eur.corr["Corr_mean"]
corr.results$Measure_Value[4] = eur2eur.corr["Corr_SD"]
corr.results$Measure_Value[5] = afr2eur.corr["Corr_mean"]
corr.results$Measure_Value[6] = afr2eur.corr["Corr_SS"]
corr.results$Measure_Value[7] = afr2afr.summary["Correlation"]
corr.results$Measure_Value[8] = 0

# compile and save to file
cat("writing results to file...\n")
fwrite(x = data.table(rbind(r2.results, corr.results)), file = file.path(eur89.resultsdir, "geuvadis.subsampling.results.bypop.summary.txt"), sep = "\t", quote = FALSE, na = "NA")

cat("\nanalysis complete.\n")
