#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -l arch=linux-x64               # -- SGE resources (CPU type)

# user limits: -c max size of core files created
date
hostname
ulimit -c 0 

genelist=$genelist
#sagegenopfx=$sagegenopfx
subjectids_eur=$subjectids_eur
subjectids_yri=$subjectids_yri
sage_lab2nwd=$sage_lab2nwd
Rscript=$Rscript
R_compute_new_weights=$R_compute_new_weights
exprfile_eur=$exprfile_eur
exprfile_yri=$exprfile_yri
logdir=$logdir
alpha=$alpha
gctadir=$gctadir
glmmethod=$glmmethod
outdir=$outdir
imputegenodir=$imputegenodir
maf=$maf
hwe=$hwe
nthreads=$nthreads
PLINK=$PLINK
memory_limit_mb=$memory_limit_mb
resultsdir=$resultsdir
resultssubdir=$resultssubdir
nwdids=$nwdids

# parse current gene
# NOTA BENE: in general BASH arrays are 0-indexed while SGE tasks are 1-indexed
# since $genelist lacks a header the genes are essentially 1-indexed
# must subtract 1 from $SGE_TASK_ID to match correct gene 
# in general, be mindful when indexing BASH arrays with SGE task IDs 
i=$(expr ${SGE_TASK_ID} - 1)
read -a genes <<< $(cat $genelist | cut -d " " -f 1)
gene=${genes[$i]}

# path to gene folder
#genedir="${gctadir}/${gene}"
#genepath="${genedir}/${gene}"
#genepath="/netapp/home/kkeys/gala_sage/rnaseq/imputed_genotypes/${gene}/${gene}"
genepath="${outdir}/${gene}"
eurgenopfx="${genepath}/${gene}_eur"
yrigenopfx="${genepath}/${gene}_eur"
mkdir -p $genepath # make gene directory in case it doesn't exist
rawpath_eur="${eurgenopfx}.raw"
rawpath_yri="${yrigenopfx}.raw"
bimfile_eur="${eurgenopfx}.bim"
bimfile_yri="${yrigenopfx}.bim"
#bimfile="/netapp/home/kkeys/gala_sage/genotypes/SAGE/mergedLAT-LATP/SAGE_mergedLAT-LATP_030816_rsID.bim"
resultssubdir="${resultsdir}/results"

# output files
#predictionfile="${outdir}/sage_predictions_${glmmethod}_${gene}.txt"
#lambdafile="${outdir}/sage_lambdas_${glmmethod}_${gene}.txt"
#weightsfile="${outdir}/sage_weights_${glmmethod}_${gene}.txt"
predictionfile="${resultssubdir}/geuvadis_predictions_${glmmethod}_${gene}.txt"
lambdafile="${resultssubdir}/geuvadis_lambdas_${glmmethod}_${gene}.txt"
weightsfile="${resultssubdir}/geuvadis_weights_${glmmethod}_${gene}.txt"

#echo "weightsfile = $weightsfile"
#echo "predictionfile = $predictionfile"
#echo "lambdafile = $lambdafile"


mygeneinfo=$(grep $gene ${genelist})
chr=$(echo $mygeneinfo | cut -f 2 -d " ")
startpos=$(echo $mygeneinfo | cut -f 3 -d " ")
endpos=$(echo $mygeneinfo | cut -f 4 -d " ")

# with $chr we point to the correct BED/BIM/BAM files
#bedfile="${imputegenodir}/SAGE_IMPUTE_HRC_chr${chr}.dose"
bedfile="${imputegenodir}/GEUVADIS.ALLCHR.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.rsq_0.8_maf_0.01_hwe_0.00001_geno_0.05"

# create a PLINK RAW file
# this codes the dosage format required for glmnet
# here we use the genome-wide GEUVADIS genotype data with rsIDs
# separate EUR and YRI into two file groups 
$PLINK \
    --bfile $bedfile \
    --chr ${chr} \
    --from-bp ${startpos} \
    --to-bp ${endpos} \
    --maf ${maf} \
    --hwe ${hwe} \
    --recode A \
    --make-bed \
    --out ${yrigenopfx} \
    --threads ${nthreads} \
    --memory ${memory_limit_mb} \
    --keep ${subjectids_yri} #\
    #--silent

$PLINK \
    --bfile $bedfile \
    --chr ${chr} \
    --from-bp ${startpos} \
    --to-bp ${endpos} \
    --maf ${maf} \
    --hwe ${hwe} \
    --recode A \
    --make-bed \
    --out ${eurgenopfx} \
    --threads ${nthreads} \
    --memory ${memory_limit_mb} \
    --keep ${subjectids_eur} #\
    #--silent

## using WGS or imputed data? the IDs probably need replacing
## replace the Subject IDs with NWDIDs one-by-one using $sage_lab2nwd
## this loop **explicitly clobbers the PLINK RAW file**
## do not make a habit of this! 
#while read -r -a line;
#do
#    nwdid="${line[0]}"
#    subjid="${line[1]}"
#    sed -i -e "s/$subjid/$nwdid/g" ${rawpath}
#done < ${sage_lab2nwd}
#
##bimfile="${sagegenopfx}.bim"
 
# call glmnet script
# the method used depends on the alpha value: 
# alpha="0.5" --> elastic net regression
# alpha="1.0" --> LASSO regression
# alpha="0.0" --> ridge regression
echo "starting R script to compute new GTEx weights..."
 $Rscript $R_compute_new_weights ${rawpath_eur} ${exprfile_eur} ${gene} ${predictionfile} ${lambdafile} ${alpha} ${weightsfile} ${bimfile_eur}
#head -n 5 ${rawpath_yri} 


#    # call R2 calculation script with default weights
#    #if [[ "$r2_default" != "0" -a "$r2_default" != ""]] ; then
#    if [[ ! -z "$r2_default" ]]; then
#        if [[ "$r2_default" != "0" ]]; then 
#            $Rscript $R_compute_r2 ${rawpath} ${exprfile} ${gene} ${ucsc_snpfile} ${r2resultsfile} "genotype_data"
#        fi
#    fi

# query return value of previous command
RETVAL=$?

# if return value is not 0, then previous command did not exit correctly
# create a status file notifying of error
# in contrary case, notify of success
if [ "$RETVAL" -ne "0" ];
then
  echo "ERROR" > ${logdir}/status.${glmmethod}.${gene}
else
  echo "SUCCESS" > ${logdir}/status.${glmmethod}.${gene}
fi

# DO NOT clean up scratch space, e.g., with `rm -rf ${scratchdir}`
# will subsequently compile results after this array job completes

#cat $predictionfile
#cat $lambdafile
#cat $weightsfile


# ask job to report on itself
echo "job ${JOB_ID} report:"
qstat -j ${JOB_ID}
