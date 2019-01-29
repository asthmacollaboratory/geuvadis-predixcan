#!/usr/bin/env bash

# genes with
# -- measurements
# -- predictions in all 5 test pops
# -- predictions across all train pops
# number: 10,631
cat geuvadis.crosspop.results.*.commongenes.txt | sort | uniq -c | grep "5 E" | wc -l


# from this reduced gene list, get imputation results and put them in a new table
# first remove previous table
# then add header
# then add data
allgenes="geuvadis.crosspop.results.allpop.to.allpop.results.txt"
allgenes_summary="geuvadis.crosspop.results.allpop.to.allpop.commongenes.summary.txt"
rm -f ${allgenes}
touch ${allgenes}
echo -e "Gene\tTest_Pop\tCorrelation_Mean\tCorrelation_pval\tR2_Mean\tR2_pval\tTrain_Pop" > ${allgenes}
grep --no-filename -F -f <(cat geuvadis.crosspop.results.*.commongenes.txt | sort | uniq -c | grep "5 E" | awk '{ print $2 }') geuvadis.crosspop.results.*.predictinto.allpop.results.txt >> ${allgenes} 

# summarize R2s
Rscript -e "library(data.table); library(dplyr); allgenes = fread(\"${allgenes}\"); allgenes %>% group_by(Train_Pop, Test_Pop) %>% summarize(R2_mean = mean(R2_Mean, na.rm = T), R2_StdErr = sd(R2_Mean, na.rm = T)) %>% as.data.table %>% fwrite(x = ., file = \"${allgenes_summary}\", sep = \"\\t\", na = 'NA', quote = F)"

# genes with
# -- measurements
# -- predictions in all 5 test pops
# -- predictions across all train pops 
# -- positive correlations for ALL train-test cases
# number: 142
cat geuvadis.crosspop.results.*.commongenes.poscorr.txt | sort | uniq -c | grep "5 E" | awk '{ print $2 }' | wc -l

# from this (tiny) gene list, get the imputation results
# stuff them into a new table
allgenes_poscorr="geuvadis.crosspop.results.allpop.to.allpop.results.poscorr.txt"
allgenes_poscorr_summary="geuvadis.crosspop.results.allpop.to.allpop.commongenes.summary.poscorr.txt"
rm -f ${allgenes_poscorr}
touch ${allgenes_poscorr}
echo -e "Gene\tTest_Pop\tCorrelation_Mean\tCorrelation_pval\tR2_Mean\tR2_pval\tTrain_Pop" > ${allgenes_poscorr}
grep --no-filename -F -f <(cat geuvadis.crosspop.results.*.commongenes.poscorr.txt | sort | uniq -c | grep "5 E" | awk '{ print $2 }') geuvadis.crosspop.results.*.predictinto.allpop.results.txt >> ${allgenes_poscorr} 

Rscript -e "library(data.table); library(dplyr); allgenes = fread(\"${allgenes_poscorr}\"); allgenes %>% group_by(Train_Pop, Test_Pop) %>% summarize(R2_mean = mean(R2_Mean, na.rm = T), R2_StdErr = sd(R2_Mean, na.rm = T)) %>% as.data.table %>% fwrite(x = ., file = \"${allgenes_poscorr_summary}\", sep = \"\\t\", na = 'NA', quote = F)"
