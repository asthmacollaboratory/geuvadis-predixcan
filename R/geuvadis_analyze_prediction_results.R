#!/usr/bin/env Rscript --vanilla
# ==============================================================================
# copyright Asthma Collaboratory (2018)
# coded by Kevin L. Keys
#
# Compile prediction results from the GEUVADIS EUR and YRI datasets.
# For each of the two populations, the results include
# -- estimated narrow-sense cis-heritabilities for each gene (h2)
# -- predictions from pop to itself
# -- prediction from pop into other pop
#
# Args: (none)
#
# Outputs:
#     * a data frame containing all merged gene-by-gene results 
#     * a data frame containing gene-by sample results for EUR
#     * a data frame containing gene-by sample results for YRI 
# ==============================================================================

# load libraries
library(data.table)
library(methods)

# script directories
home       = Sys.getenv("HOME")
resultsdir = paste(home, "/gala_sage/rnaseq/glmnet/elasticnet", sep = "")
rna.dir    = paste(home, "/gala_sage/rnaseq/geuvadis", sep = "")
eur.dir    = paste(resultsdir, "/eur373", sep = "")
yri.dir    = paste(resultsdir, "/yri89",  sep = "")

# script file paths
eur.pred.file    = paste(eur.dir, "/geuvadis_elasticnet_eur373_predictions.txt", sep = "")
yri.pred.file    = paste(yri.dir, "/geuvadis_elasticnet_yri89_predictions.txt",  sep = "")
eur.meas.file    = paste(rna.dir, "/geuvadis.eur373.RPKM.invnorm.txt", sep = "") 
yri.meas.file    = paste(rna.dir, "/geuvadis.yri89.RPKM.invnorm.txt",  sep = "") 
eur.h2.file      = paste(eur.dir, "/geuvadis_h2_eur373.txt", sep = "")
yri.h2.file      = paste(yri.dir, "/geuvadis_h2_yri89.txt",  sep = "")
eur.numpred.file = paste(eur.dir, "/geuvadis_elasticnet_eur373_numpred.txt", sep = "")
yri.numpred.file = paste(yri.dir, "/geuvadis_elasticnet_yri89_numpred.txt", sep = "")
eur.merge.file   = paste(resultsdir, "/geuvadis_eur373_results.txt", sep = "") 
yri.merge.file   = paste(resultsdir, "/geuvadis_yri89_results.txt",  sep = "") 
eur.h2.null.file = paste(eur.dir, "/geuvadis_h2_null_eur373.txt", sep = "")
yri.h2.null.file = paste(yri.dir, "/geuvadis_h2_null_yri89.txt",  sep = "")
eur.genelm.file  = paste(eur.dir, "/geuvadis_elasticnet_eur373_genelm_predvmeas_results.txt", sep = "")
yri.genelm.file  = paste(yri.dir, "/geuvadis_elasticnet_yri89_genelm_predvmeas_results.txt",  sep = "")

genelm.merge.file        = paste(resultsdir, "/geuvadis_genelm_results.txt",  sep = "") 
yri.pred.from.eur.file   = paste(eur.dir, "/geuvadis_elasticnet_eur373_predictinto_yri89.txt", sep = "")
eur.pred.from.yri.file   = paste(yri.dir, "/geuvadis_elasticnet_yri89_predictinto_eur373.txt", sep = "")
yri.from.eur.genelm.file = paste(eur.dir, "/geuvadis_elasticnet_eur373_predictinto_yri89_genelm_predvmeas_results.txt", sep = "") 
eur.from.yri.genelm.file = paste(yri.dir, "/geuvadis_elasticnet_yri89_predictinto_eur373_genelm_predvmeas_results.txt", sep = "") 

# load data
# load gene-subject data first
eur.pred = fread(eur.pred.file)
eur.meas = fread(eur.meas.file)
yri.pred = fread(yri.pred.file)
yri.meas = fread(yri.meas.file)

eur.pred.from.yri = fread(eur.pred.from.yri.file)
yri.pred.from.eur = fread(yri.pred.from.eur.file)

# load gene-based data next
eur.h2      = fread(eur.h2.file)
eur.null.h2 = fread(eur.h2.file)
eur.genelm  = fread(eur.genelm.file)
eur.numpred = fread(eur.numpred.file)
yri.h2      = fread(yri.h2.file)
yri.null.h2 = fread(yri.h2.file)
yri.genelm  = fread(yri.genelm.file)
yri.numpred = fread(yri.numpred.file)

eur.from.yri.genelm = fread(eur.from.yri.genelm.file, na.strings = c("NA", "NA "))
yri.from.eur.genelm = fread(yri.from.eur.genelm.file, na.strings = c("NA", "NA "))

# need header for yri.pred.from.eur 
my.header.yri = c("Gene","NA18486","NA18487","NA18488","NA18489","NA18498","NA18499","NA18502","NA18505","NA18508","NA18510","NA18511","NA18517","NA18519","NA18520","NA18858","NA18861","NA18867","NA18868","NA18870","NA18873","NA18907","NA18908","NA18909","NA18910","NA18912","NA18916","NA18917","NA18923", "NA18933","NA18934","NA19092","NA19093","NA19095","NA19096","NA19098","NA19099","NA19102","NA19107","NA19108","NA19113","NA19114","NA19116","NA19117","NA19118","NA19119","NA19121","NA19129","NA19130","NA19131","NA19137","NA19138","NA19141","NA19143","NA19144","NA19146","NA19147","NA19149", "NA19150","NA19152","NA19153","NA19159","NA19160","NA19171","NA19172","NA19175","NA19184","NA19185","NA19189","NA19190","NA19197","NA19198","NA19200","NA19201","NA19204","NA19206","NA19207","NA19209","NA19210","NA19213","NA19214","NA19222","NA19223","NA19225","NA19235","NA19236","NA19247","NA19248","NA19256","NA19257")
names(yri.pred) = my.header.yri

# do same for eur.pred.from.yri
my.header.eur =
c("Gene","HG00096","HG00097","HG00099","HG00100","HG00101","HG00102","HG00103","HG00104","HG00105","HG00106","HG00108","HG00109","HG00110","HG00111","HG00112","HG00114","HG00115","HG00116","HG00117","HG00118","HG00119","HG00120","HG00121","HG00122","HG00123","HG00124","HG00125","HG00126","HG00127","HG00128","HG00129","HG00130","HG00131","HG00132","HG00133","HG00134","HG00135","HG00136","HG00137","HG00138","HG00139","HG00141","HG00142","HG00143","HG00145","HG00146","HG00148","HG00149","HG00150","HG00151","HG00152","HG00154","HG00155","HG00156","HG00157","HG00158","HG00159","HG00160","HG00171","HG00173","HG00174","HG00176","HG00177","HG00178","HG00179","HG00180","HG00181","HG00182","HG00183","HG00185","HG00186","HG00187","HG00188","HG00189","HG00231","HG00232","HG00233","HG00234","HG00235","HG00236","HG00238","HG00239","HG00240","HG00242","HG00243","HG00244","HG00245","HG00246","HG00247","HG00249","HG00250","HG00251","HG00252","HG00253","HG00255","HG00256","HG00257","HG00258","HG00259","HG00260","HG00261","HG00262","HG00263","HG00264","HG00265","HG00266","HG00267","HG00268","HG00269","HG00271","HG00272","HG00273","HG00274","HG00275","HG00276","HG00277","HG00278","HG00280","HG00281","HG00282","HG00284","HG00285","HG00306","HG00308","HG00309","HG00310","HG00311","HG00312","HG00313","HG00315","HG00319","HG00320","HG00321","HG00323","HG00324","HG00325","HG00326","HG00327","HG00328","HG00329","HG00330","HG00331","HG00332","HG00334","HG00335","HG00336","HG00337","HG00338","HG00339","HG00341","HG00342","HG00343","HG00344","HG00345","HG00346","HG00349","HG00350","HG00351","HG00353","HG00355","HG00356","HG00358","HG00359","HG00360","HG00361","HG00362","HG00364","HG00365","HG00366","HG00367","HG00369","HG00371","HG00372","HG00373","HG00375","HG00376","HG00377","HG00378","HG00379","HG00380","HG00381","HG00382","HG00383","HG00384","HG01334","HG01789","HG01790","HG01791","HG02215","NA06984","NA06985","NA06986","NA06989","NA06994","NA07037","NA07048","NA07051","NA07056","NA07346","NA07347","NA07357","NA10847","NA10851","NA11829","NA11830","NA11831","NA11832","NA11840","NA11843","NA11881","NA11892","NA11893","NA11894","NA11918","NA11920","NA11930","NA11931","NA11992","NA11993","NA11994","NA11995","NA12004","NA12005","NA12006","NA12043","NA12044","NA12045","NA12058","NA12144","NA12154","NA12155","NA12156","NA12234","NA12249","NA12272","NA12273","NA12275","NA12282","NA12283","NA12286","NA12287","NA12340","NA12341","NA12342","NA12347","NA12348","NA12383","NA12399","NA12400","NA12413","NA12489","NA12546","NA12716","NA12717","NA12718","NA12749","NA12750","NA12751","NA12760","NA12761","NA12762","NA12763","NA12775","NA12776","NA12777","NA12778","NA12812","NA12813","NA12814","NA12815","NA12827","NA12829","NA12830","NA12842","NA12843","NA12872","NA12873","NA12874","NA12889","NA12890","NA20502","NA20503","NA20504","NA20505","NA20506","NA20507","NA20508","NA20509","NA20510","NA20512","NA20513","NA20514","NA20515","NA20516","NA20517","NA20518","NA20519","NA20520","NA20521","NA20524","NA20525","NA20527","NA20528","NA20529","NA20530","NA20531","NA20532","NA20534","NA20535","NA20536","NA20537","NA20538","NA20539","NA20540","NA20541","NA20542","NA20543","NA20544","NA20581","NA20582","NA20585","NA20586","NA20588","NA20589","NA20752","NA20754","NA20756","NA20757","NA20758","NA20759","NA20760","NA20761","NA20765","NA20766","NA20768","NA20769","NA20770","NA20771","NA20772","NA20773","NA20774","NA20778","NA20783","NA20785","NA20786","NA20787","NA20790","NA20792","NA20795","NA20796","NA20797","NA20798","NA20799","NA20800","NA20801","NA20802","NA20803","NA20804","NA20805","NA20806","NA20807","NA20808","NA20809","NA20810","NA20811","NA20812","NA20813","NA20814","NA20815","NA20816","NA20819","NA20826","NA20828") 
names(eur.pred) = my.header.eur

# melt data tables and rename column names for YRI 
yri.meas.melt = melt(yri.meas, id.var = "Gene")
yri.pred.melt = melt(yri.pred, id.var = "Gene")
colnames(yri.meas.melt) = c("Gene", "SubjectID", "Measured_Expr")
colnames(yri.pred.melt) = c("Gene", "SubjectID", "Predicted_Expr")
colnames(yri.h2)        = c("Gene", "h2_YRI", "h2_stderr_YRI")
colnames(yri.null.h2)   = c("Gene", "h2_null_YRI", "h2_null_stderr_YRI")
colnames(yri.genelm)    = c("Gene", "P.value_YRI", "R2_YRI", "Spearman.rho_YRI")
colnames(yri.numpred)   = c("Gene", "Num_Pred_YRI")

# do same for EUR
eur.pred.melt = melt(eur.pred, id.var = "Gene")
eur.meas.melt = melt(eur.meas, id.var = "Gene")
colnames(eur.meas.melt) = c("Gene", "SubjectID", "Measured_Expr")
colnames(eur.pred.melt) = c("Gene", "SubjectID", "Predicted_Expr")
colnames(eur.h2)        = c("Gene", "h2_EUR", "h2_stderr_EUR")
colnames(eur.null.h2)   = c("Gene", "h2_null_EUR", "h2_null_stderr_EUR")
colnames(eur.genelm)    = c("Gene", "P.value_EUR", "R2_EUR", "Spearman.rho_EUR")
colnames(eur.numpred)   = c("Gene", "Num_Pred_EUR")

# handle crosspopulation predictions
colnames(yri.pred.from.eur)   = c("Gene", "SubjectID", "Predicted_Expr_YRIfromEUR")
colnames(eur.pred.from.yri)   = c("Gene", "SubjectID", "Predicted_Expr_EURfromYRI")
colnames(yri.from.eur.genelm) = c("Gene", "P.value_YRIfromEUR", "R2_YRIfromEUR", "Spearman_rho_YRIfromEUR")
colnames(eur.from.yri.genelm) = c("Gene", "P.value_EURfromYRI", "R2_EURfromYRI", "Spearman_rho_EURfromYRI")

# merge YRI gene-subject results 
yri.merge = merge(yri.meas.melt, yri.pred.melt, by = c("Gene", "SubjectID"), all = TRUE)
yri.merge = merge(yri.merge, yri.pred.from.eur, by = c("Gene", "SubjectID"), all.x = TRUE)

# merge EUR gene-subject results
eur.merge = merge(eur.meas.melt, eur.pred.melt, by = c("Gene", "SubjectID"), all = TRUE) 
eur.merge = merge(eur.merge, eur.pred.from.yri, by = c("Gene", "SubjectID"), all = TRUE) 

# now merge genewise results
genelm.merge = merge(eur.genelm,   yri.genelm,  by = c("Gene"), all = TRUE)
genelm.merge = merge(genelm.merge, eur.h2,      by = c("Gene"), all = TRUE)
genelm.merge = merge(genelm.merge, yri.h2,      by = c("Gene"), all = TRUE)
genelm.merge = merge(genelm.merge, eur.null.h2, by = c("Gene"), all = TRUE)
genelm.merge = merge(genelm.merge, yri.null.h2, by = c("Gene"), all = TRUE)
genelm.merge = merge(genelm.merge, eur.numpred, by = c("Gene"), all.x = TRUE)
genelm.merge = merge(genelm.merge, yri.numpred, by = c("Gene"), all.x = TRUE)
genelm.merge = merge(genelm.merge, yri.from.eur.genelm, by = c("Gene"), all.x = TRUE)
genelm.merge = merge(genelm.merge, eur.from.yri.genelm, by = c("Gene"), all.x = TRUE)

# save merged files
fwrite(x = eur.merge, file = eur.merge.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", na = "NA")
fwrite(x = yri.merge, file = yri.merge.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", na = "NA")
fwrite(x = genelm.merge, file = genelm.merge.file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t", na = "NA")
