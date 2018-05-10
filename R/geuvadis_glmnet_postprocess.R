#!/usr/bin/env Rscript --vanilla

suppressMessages(library(dplyr))
suppressMessages(library(data.table))


#geuvadis.glmnet.postprocess = function(){

    # read command line arguments
    args = commandArgs(trailingOnly = TRUE)

    # parse command line arguments
    weightsfile            = args[1]
    new.weightsfile        = args[2]
    discard.ratio          = as.numeric(args[3])
    num.pred.file          = args[4]
    num.samples            = as.integer(args[5])
    prediction.file        = args[6]
    exprfile               = args[7]
    out.lm.file            = args[8]
    out.genelm.file        = args[9]
    altpop.pred.file       = args[10]
    altpop.exprfile        = args[11]
    altpop.out.lm.file     = args[12]
    altpop.out.genelm.file = args[13]

    # open weights file
    x = fread(weightsfile)

    # ensure proper header
    my.header = c("Gene","Held_out_Sample","SNP","A1","A2","Beta")
    colnames(x) = my.header

    # subset for weights with no NA values
#    x.nona = x[!is.na(x$Beta),]

    # determine how many predicted samples that each gene has
    # we will sink this into an output file and reload it as a data.table
    cat(paste("creating num.pred.file at ", num.pred.file, "...\n", sep = ""))
#    sink(num.pred.file)
#    for (i in unique(sort(x$Gene))){
#
#        # output format is "(gene)\t(#pred samples)"
#        #cat(i, "\t", length(unique(sort(x.nona[x.nona$Gene == i,]$Held_out_Sample))), "\n")
#        cat(i, "\t", length(unique(sort(x[x$Gene == i,]$Held_out_Sample))), "\n")
#    }
#    sink()
    num.pred = data.table(count(count(x, Gene, Held_out_Sample), Gene))
    colnames(num.pred) = c("Gene", "Num.Pred")
    fwrite(x = num.pred, file = num.pred.file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t") 
    cat("done.\n")

    # reload the formatted file
    #cat("now reading num.pred.file from disk\n")
    #num.pred = fread(num.pred.file)

    # subset those genes with at least discard.ratio * (# samples) predicted samples
    cat("discarding genes with less than ", discard.ratio * num.samples, "predictions...\n")
    num.pred.sub = num.pred[num.pred$Num.Pred > discard.ratio*num.samples,]

    # we can now look at the predictions for the "well-predicted" genes in num.pred.sub
    # must load prediction information first
    cat("reading prediction.file at ", prediction.file, "...")
    geuvadis.predictions = fread(prediction.file)

    # now load measured RNA data
    cat("reading gene expression file at ", exprfile, "...\n")
    geuvadis.rnaseq = fread(exprfile)

    # expected output, using YRI 89:
    # > colnames(geuvadis.rnaseq)
    #     [1] "Gene"    "NA18486" "NA18487" "NA18488" "NA18489" "NA18498" "NA18499"
    #     [8] "NA18502" "NA18505" "NA18508" "NA18510" "NA18511" "NA18517" "NA18519"
    #    [15] "NA18520" "NA18858" "NA18861" "NA18867" "NA18868" "NA18870" "NA18873"
    #    [22] "NA18907" "NA18908" "NA18909" "NA18910" "NA18912" "NA18916" "NA18917"
    #    [29] "NA18923" "NA18933" "NA18934" "NA19092" "NA19093" "NA19095" "NA19096"
    #    [36] "NA19098" "NA19099" "NA19102" "NA19107" "NA19108" "NA19113" "NA19114"
    #    [43] "NA19116" "NA19117" "NA19118" "NA19119" "NA19121" "NA19129" "NA19130"
    #    [50] "NA19131" "NA19137" "NA19138" "NA19141" "NA19143" "NA19144" "NA19146"
    #    [57] "NA19147" "NA19149" "NA19150" "NA19152" "NA19153" "NA19159" "NA19160"
    #    [64] "NA19171" "NA19172" "NA19175" "NA19184" "NA19185" "NA19189" "NA19190"
    #    [71] "NA19197" "NA19198" "NA19200" "NA19201" "NA19204" "NA19206" "NA19207"
    #    [78] "NA19209" "NA19210" "NA19213" "NA19214" "NA19222" "NA19223" "NA19225"
    #    [85] "NA19235" "NA19236" "NA19247" "NA19248" "NA19256" "NA19257"

    # the rnaseq and prediction files have the same column name order
    colnames(geuvadis.predictions) = colnames(geuvadis.rnaseq)

    # now subset predictions based on which genes are well predicted
    #colnames(geuvadis.predictions)[1] = "Gene"
    geuvadis.predictions.sub = geuvadis.predictions[geuvadis.predictions$Gene %in% num.pred.sub$Gene,]


    # subset the file with just the genes from geuvadis.predictions.sub
    cat("subsetting expression file to well-predicted genes...\n") 
    geuvadis.rnaseq.sub = geuvadis.rnaseq[geuvadis.rnaseq$Gene %in% geuvadis.predictions.sub$Gene,]

    # sort the genes before transposing
    cat("sorting and transposing expressions...\n")
    setorder(geuvadis.rnaseq.sub, Gene)
    #trnaseq = data.table(t(geuvadis.rnaseq.sub[,-c(1:4)]))
    trnaseq = data.table(t(geuvadis.rnaseq.sub[,-c(1)]))
    colnames(trnaseq) = geuvadis.rnaseq.sub$Gene
    trnaseq = cbind(colnames(geuvadis.rnaseq.sub)[-1], trnaseq)
    colnames(trnaseq)[1] = "SubjectID"

    # put column names and a column of sample names
    # make sure to sort by sample
    #trnaseq = cbind(my.header[-c(1:4)], trnaseq)
    setorder(trnaseq, SubjectID)

    # melt the predicted and measured expression data.tables
    cat("melting both prediction and measured expression data tables...\n")
    pred.melted   = melt(geuvadis.predictions.sub, id.vars = c("Gene"))
    #rnaseq.melted = melt(geuvadis.rnaseq.sub[,-c(1:3)], id.vars = "Gene")
    rnaseq.melted = melt(geuvadis.rnaseq.sub, id.vars = "Gene")

    # standardize their column names too
    colnames(pred.melted)   = c("Gene", "SubjectID", "Predicted_Expr")
    colnames(rnaseq.melted) = c("Gene", "SubjectID", "Measured_Expr")

    # now merge the data frames
    cat("merging melted data tables...")
    geuvadis.rnapred = merge(pred.melted, rnaseq.melted, by = c("Gene", "SubjectID"))
    setorder(geuvadis.rnapred, Gene, SubjectID)

    # do two linear models
    # first regress all (Gene, SubjectID) pairs of predicted expression onto measured expression
    cat("first linear model: regress predicted onto measured expression for all (gene, subject) pairs\n") 
    sink(out.lm.file)
    print(summary(lm(Predicted_Expr ~ Measured_Expr, data = geuvadis.rnapred)))
    sink()

    # now do individual gene regressions
    cat("second linear model: individual regressions for each gene\n")
    sink(out.genelm.file)
    cat("Gene\tP.value\tR2\tSpearman.rho\n")
    for(gene in unique(sort(geuvadis.rnapred$Gene))){
        geuvadis.rnapred.sub = geuvadis.rnapred[geuvadis.rnapred$Gene == gene,]
        my.lm = summary(lm(Predicted_Expr ~ Measured_Expr, data = geuvadis.rnapred.sub))
        lm.p = my.lm$coefficients[2,4]
        lm.r2 = my.lm$r.squared
        my.rho = cor.test(geuvadis.rnapred.sub$Predicted_Expr, geuvadis.rnapred.sub$Measured_Expr, method = "spearman")$estimate
        cat(gene, "\t", lm.p, "\t", lm.r2, "\t", my.rho, "\n")
    }
    sink()

    # now let us compare crosspopulation predictions 
    cat("now analyzing crosspopulation predictions")
    altpop.predictions = fread(altpop.pred.file)
    altpop.expr        = fread(altpop.exprfile)

    # melt the data.table for the measured expression data
    altpop.expr.melt = melt(data = altpop.expr, id.vars = c("Gene"))
    colnames(altpop.expr.melt) = c("Gene", "SubjectID", "Measured_Expr")

    # rename column names for altpop.predictions
    colnames(altpop.predictions) = c("Gene", "SubjectID", "Predicted_Expr")

    # merge the measured and predicted expressions in the alternative population 
    altpop.rnapred = merge(altpop.expr.melt, altpop.predictions, by = c("Gene", "SubjectID"))

    # first model: regress measured, predicted expression in all (gene, subject) pairs
    cat("first linear model in other population: regress predicted onto measured expression for all (gene, subject) pairs\n")  
    sink(altpop.out.lm.file)
    print(summary(lm(Predicted_Expr ~ Measured_Expr, data = altpop.rnapred)))
    sink()

    cat("second linear model in other population: individual regressions for each gene\n")
    sink(altpop.out.genelm.file)
    cat("Gene\tP.value\tR2\tSpearman.rho\n")
    for(gene in unique(sort(altpop.rnapred$Gene))){
        geuvadis.rnapred.sub = altpop.rnapred[altpop.rnapred$Gene == gene,]
        if(all(is.na(geuvadis.rnapred.sub$Predicted_Expr))){
            cat(gene, "\tNA\tNA\tNA\n")
        } else {
        my.lm = summary(lm(Predicted_Expr ~ Measured_Expr, data = geuvadis.rnapred.sub))
        lm.p = my.lm$coefficients[2,4]
        lm.r2 = my.lm$r.squared
        my.rho = cor.test(geuvadis.rnapred.sub$Predicted_Expr, geuvadis.rnapred.sub$Measured_Expr, method = "spearman")$estimate
        cat(gene, "\t", lm.p, "\t", lm.r2, "\t", my.rho, "\n")
        }
    }
    sink()


    cat("all done!\n\n")

#    return()
#}
#
#geuvadis.glmnet.postprocess()

# show any warnings
cat("any warnings?\n")
warnings()

