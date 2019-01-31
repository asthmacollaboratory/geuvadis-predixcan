#!/usr/bin/env Rscript --vanilla
#
# coded by Kevin L. Keys (2018)

# load libraries
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(dunn.test))
suppressMessages(library(methods))
suppressMessages(library(ggplot2))

# plotting options
options(bitmapType='cairo')


data.dir = "~/gala_sage/rnaseq/glmnet/elasticnet/crosspop/crosspop_results/2019-01-25"

# format Utahn data
ceu = fread(file.path(data.dir, "geuvadis.crosspop.results.ceu.predictinto.allpop.results.txt"), header = TRUE)
ceu.poscorr.transcripts = fread(file.path(data.dir, "geuvadis.crosspop.results.ceu.commongenes.poscor.txt"), header = FALSE)[[1]]
colnames(ceu)[c(2,3,5)] = c("Test_Pop", "R2", "Correlation")
ceu$Train_Pop = "CEU"
ceu.sub = ceu %>% na.omit %>% dplyr::filter(Gene %in% ceu.poscorr.transcripts) %>% select(Train_Pop, Test_Pop, Gene, R2, Correlation) %>% as.data.table

# format British data
gbr = fread(file.path(data.dir, "geuvadis.crosspop.results.gbr.predictinto.allpop.results.txt"), header = TRUE)
gbr.poscorr.transcripts = fread(file.path(data.dir, "geuvadis.crosspop.results.gbr.commongenes.poscor.txt"), header = FALSE)[[1]]
colnames(gbr)[c(2,3,5)] = c("Test_Pop", "R2", "Correlation")
gbr$Train_Pop = "GBR"
gbr.sub = gbr %>% na.omit %>% dplyr::filter(Gene %in% gbr.poscorr.transcripts) %>% select(Train_Pop, Test_Pop, Gene, R2, Correlation) %>% as.data.table

# format Tuscan data
tsi = fread(file.path(data.dir, "geuvadis.crosspop.results.tsi.predictinto.allpop.results.txt"), header = TRUE)
tsi.poscorr.transcripts = fread(file.path(data.dir, "geuvadis.crosspop.results.tsi.commongenes.poscor.txt"), header = FALSE)[[1]]
tsi$Train_Pop = "TSI"
colnames(tsi)[c(2,3,5)] = c("Test_Pop", "R2", "Correlation")
tsi.sub = tsi %>% na.omit %>% dplyr::filter(Gene %in% tsi.poscorr.transcripts) %>% select(Train_Pop, Test_Pop, Gene, R2, Correlation) %>% as.data.table

# format Finnish data
fin = fread(file.path(data.dir, "geuvadis.crosspop.results.fin.predictinto.allpop.results.txt"), header = TRUE)
fin.poscorr.transcripts = fread(file.path(data.dir, "geuvadis.crosspop.results.fin.commongenes.poscor.txt"), header = FALSE)[[1]]
fin$Train_Pop = "FIN"
colnames(fin)[c(2,3,5)] = c("Test_Pop", "R2", "Correlation")
fin.sub = fin %>% na.omit %>% dplyr::filter(Gene %in% fin.poscorr.transcripts) %>% select(Train_Pop, Test_Pop, Gene, R2, Correlation) %>% as.data.table

# format Yoruba data
yri = fread(file.path(data.dir, "geuvadis.crosspop.results.yri.predictinto.allpop.results.txt"), header = TRUE)
yri.poscorr.transcripts = fread(file.path(data.dir,"geuvadis.crosspop.results.yri.commongenes.poscor.txt"), header = FALSE)[[1]]
yri$Train_Pop = "YRI"
colnames(yri)[c(2,3,5)] = c("Test_Pop", "R2", "Correlation")
yri.sub = yri %>% na.omit %>% dplyr::filter(Gene %in% yri.poscorr.transcripts) %>% select(Train_Pop, Test_Pop, Gene, R2, Correlation) %>% as.data.table


# combine data into single data table
allpop = rbindlist(list(ceu.sub, gbr.sub, tsi.sub, fin.sub, yri.sub))

# add labels to distinguish various crosspop scenarios:
# EUR_AFR = EUR if train-test pops are both EUR, AFR otherwise
# Train_Test includes the *training* population, followed by the *testing* population
allpop = allpop %>%
    mutate(
        Train_Test = paste(Train_Pop, Test_Pop, sep = "_"),
        EUR_AFR = paste(
            ifelse(Train_Pop == "YRI" , "AFR", "EUR"),
            ifelse(Test_Pop == "YRI", "AFR", "EUR"),
            sep = "_"
        )
    ) %>%
    as.data.table

# this data.table houses just the cross-population imputation results
allpop.sub = allpop %>% dplyr::filter(Train_Pop != Test_Pop) %>% as.data.table 


# this data.table only has results for genes in common to all continental population testing scenarios
commongenes.eurafr = allpop %>%
    count(EUR_AFR, Gene) %>%
    select(Gene) %>%
    count(Gene) %>%
    dplyr::filter(n == 4) %>%
    select(Gene) %>%
    unlist %>%
    unname
allpop.eurafr = allpop %>% dplyr::filter(Gene %in% commongenes.eurafr) %>% as.data.table

# this data.table only has results for genes in common to all train-test scenarios
commongenes.traintest = allpop %>%
    count(Train_Test, Gene) %>%
    select(Gene) %>%
    count(Gene) %>%
    dplyr::filter(n == 25) %>%
    select(Gene) %>%
    unlist %>%
    unname
allpop.traintest = allpop %>% dplyr::filter(Gene %in% commongenes.traintest) %>% as.data.table

# Kruskal-Wallis tests determine significant differences in means 
#kruskal.eurafr = kruskal.test(R2 ~ as.factor(EUR_AFR),    data = allpop.sub)
#kruskal.test(R2 ~ as.factor(Train_Test), data = allpop.sub)
#kruskal.test(R2 ~ as.factor(Train_Pop),  data = allpop.sub)
#kruskal.test(R2 ~ as.factor(Test_Pop),   data = allpop.sub)

# Dunn tests tease apart which groups are sig. diff. by KW test
# output list contains $P.adjusted, $Z, and $comparisons fields to describe p-values and KW statistics
r2.by.eurafr.results = data.frame(dunn.test(x = allpop.eurafr$R2, g = allpop.eurafr$EUR_AFR, method = "bonferroni")[-1])
#r2.by.eurafr = dunn.test(x = allpop.sub$R2, g = allpop.sub$EUR_AFR, method = "bonferroni")
#dunn.test(x = allpop.sub$R2, g = allpop.sub$Train_Pop,  method = "bonferroni")
#dunn.test(x = allpop.sub$R2, g = allpop.sub$Test_Pop,   method = "bonferroni")
#r2.by.eurafr = dunn.test(x = allpop.sub$R2, g = allpop.sub$Train_Test, method = "bonferroni")  ## this takes awhile to compute
r2.by.traintest.results = data.frame(dunn.test(x = allpop.traintest$R2, g = allpop.traintest$Train_Test, method = "bonferroni")[-1])  ## this takes awhile to compute

# make boxplots of results using EUR_AFR
g.r2byeurafr = ggplot(allpop, aes(x = as.factor(EUR_AFR), y = R2), color = EUR_AFR) +
    geom_point() +
    geom_boxplot(alpha = 0.8) +
    ggtitle(bquote("Comparison of imputation " ~ R^2 ~ " between continental GEUVADIS populations")) +
    xlab("Train pop - Test Pop") +
    ylab(bquote(R^2))
ggsave(g.r2byeurafr, filename = "geuvadis.crosspop.results.R2.by.eurafr.png")

g.r2byeurafr.commongenes = ggplot(allpop.eurafr, aes(x = as.factor(EUR_AFR), y = R2), color = EUR_AFR) +
    geom_point() +
    geom_boxplot(alpha = 0.8) +
    ggtitle(bquote("Comparison of " ~ R^2 ~ " between continental GEUVADIS populations"), subtitle = bquote("N = " ~ .(length(commongenes.eurafr)) ~ " common genes")) +
    xlab("Train pop - Test pop") +
    ylab(bquote(R^2))
ggsave(g.r2byeurafr.commongenes, filename = "geuvadis.crosspop.results.R2.by.eurafr.commongenes.png")

# make boxplots of results using Train_Test 
g.r2byeurafr = ggplot(allpop, aes(x = as.factor(Train_Test), y = R2), color = Train_Test) +
    geom_point() +
    geom_boxplot(alpha = 0.8) +
    ggtitle(bquote("Comparison of " ~ R^2 ~ " between all GEUVADIS populations")) +
    xlab("Train pop - Test pop") +
    ylab(bquote(R^2))
ggsave(g.r2byeurafr, filename = "geuvadis.crosspop.results.R2.by.traintest.png", units = "in", width = 20, height = 7)

g.r2byeurafr.commongenes = ggplot(allpop.traintest, aes(x = as.factor(Train_Test), y = R2), color = Train_Test) +
    geom_point() +
    geom_boxplot(alpha = 0.8) +
    ggtitle(bquote("Comparison of " ~ R^2 ~ " between all GEUVADIS populations"), subtitle = bquote("N = " ~ .(length(commongenes.traintest)) ~ " common genes")) +
    xlab("Train pop - Test pop") +
    ylab(bquote(R^2))

ggsave(g.r2byeurafr.commongenes, filename = "geuvadis.crosspop.results.R2.by.traintest.commongenes.png", units = "in", width = 20, height = 7)
