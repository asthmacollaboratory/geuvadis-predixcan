#!/usr/bin/env Rscript --vanilla
library(data.table)
library(dplyr)

set.seed(2018)

# output file path names
yri2eur.nofinn.out.genelm.file = "geuvadis_elasticnet_yri89_predictinto_eur278nofinn_genelm_predvmeas_results.txt"   
yri2eur89.out.genelm.file = "geuvadis_elasticnet_yri89_predictinto_eur89_genelm_predvmeas_results.txt" 

# read files
yri2eur = fread("../yri89/geuvadis_elasticnet_yri89_predictinto_eur373.txt")
finn = fread("~/gala_sage/rnaseq/geuvadis/geuvadis.fin95.sampleids.txt", header = FALSE)
eur.rnaseq = fread("~/gala_sage/rnaseq/geuvadis/geuvadis.eur373.RPKM.invnorm.txt")

# melt and rename the EUR RNA-seq data in preparation for merging
# then merge with the YRI -> EUR predictions
# finally, rename the prediction and measurement column for ease of lm()-ing
eur.rnaseq.df = melt(eur.rnaseq)
colnames(eur.rnaseq.df) = c("Gene", "SubjectID", "Measured_Expr")
yri2eur = merge(yri2eur, eur.rnaseq.df, by = c("Gene", "SubjectID"), all = TRUE)
colnames(yri2eur) = c("Gene", "SubjectID", "Predicted_Expr", "Measured_Expr")

# make a subset of EUR with no Finnish individuals included
# create a rectangular data frame in prep for lm()-ing
yri2eur.nofinn = yri2eur[!(yri2eur$SubjectID %in% finn$V1), ]
yri2eur.nofinn.df = dcast(data = yri2eur.nofinn, formula = SubjectID ~ Gene, value.var = "Prediction_Expr")

# make a random subset of EUR that is equal in size to YRI training set
# pull those samples from the prediction results
# then make a rectangular data frame as before in prep for lm()-ing
eur89.ids = sample(levels(factor(yri2eur.nofinn$SubjectID)), size = 89, replace = FALSE)
yri2eur89 = yri2eur[(yri2eur$SubjectID %in% finn$V1), ]
yri2eur89.df = dcast(data = yri2eur89, formula = SubjectID ~ Gene, value.var = "Prediction_Expr")






compute.gene.lm = function(data, out.file, genes = unique(sort(data$Gene))){
    #sink(out.file)
    #cat("Gene\tP.value\tR2\tSpearman.rho\n")
    ngenes = length(genes)
    out.matrix = matrix(NA, ngenes, 4)
    out.matrix[,1] = genes
    #for(gene in genes){
    for (i in 1:ngenes) {
        gene = genes[i]
        data.sub = data[data$Gene == gene,]
        #if(all(is.na(data.sub$Predicted_Expr))){
            #cat(gene, "\tNA\tNA\tNA\n")
        #} else {
        if(!all(is.na(data.sub$Predicted_Expr))){
            my.lm  = summary(lm(Predicted_Expr ~ Measured_Expr, data = data.sub))
            lm.p   = my.lm$coefficients[2,4]
            lm.r2  = my.lm$r.squared
            my.rho = cor.test(data.sub$Predicted_Expr, data.sub$Measured_Expr, method = "spearman")$estimate
        #cat(gene, "\t", lm.p, "\t", lm.r2, "\t", my.rho, "\n")
            out.matrix[i,2:4] = c(lm.p, lm.r2, my.rho)
        }   
    }   
    #sink()
    out.df = data.table(out.matrix)
    colnames(out.df) = c("Gene", "P.value", "R2", "Spearman.rho")
    fwrite(out.df, file = out.file, row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA", sep = "\t")
    return()
}

#sink(yri2eur.nofinn.out.genelm.file)
#compute.gene.lm(yri2eur.nofinn)
#sink()
genes = unique(sort(yri2eur.nofinn$Gene))
compute.gene.lm(yri2eur.nofinn, yri2eur.nofinn.out.genelm.file, genes = genes)
genes = unique(sort(yri2eur89$Gene))
compute.gene.lm(yri2eur89, yri2eur89.out.genelm.file, genes = genes)
