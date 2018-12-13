#!/usr/bin/env bash                # -- what is the language of this shell?
#                                  # -- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                    # -- the shell for the job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -l arch=linux-x64               # -- SGE resources (CPU type)

# user limits: -c max size of core files created
echo "Started script at timestamp $(date)"
hostname
ulimit -c 0 

genelist=$genelist
subjectids=$subjectids
Rscript=$Rscript
R_compute_new_weights=$R_compute_new_weights
exprfile=$exprfile
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
pop=$pop
altpop=$altpop
tmpdir=$tmpdir
subjectids_altpop=$subjectids_altpop
R_predict_new_pop=$R_predict_new_pop
GCTA=$GCTA
phenofile=$phenofile
counter=$SGE_TASK_ID
FIESTA=$FIESTA
PYTHON=$PYTHON
nfolds=$nfolds


# parse current gene
# NOTA BENE: in general BASH arrays are 0-indexed while SGE tasks are 1-indexed
# since $genelist lacks a header the genes are essentially 1-indexed
# must subtract 1 from $SGE_TASK_ID to match correct gene 
# in general, be mindful when indexing BASH arrays with SGE task IDs 
i=$(expr ${SGE_TASK_ID} - 1)
read -a genes <<< $(cat $genelist | cut -d " " -f 1)
gene=${genes[$i]}

echo "Preparing analysis of ${gene} in population ${pop}..."

# path to gene folder
genepath="${outdir}/${gene}"
genopfx="${genepath}/${gene}"
mkdir -p $genepath # make gene directory in case it doesn't exist
rawpath="${genopfx}.raw"
bimfile="${genopfx}.bim"
#bimfile="/netapp/home/kkeys/gala_sage/genotypes/SAGE/mergedLAT-LATP/SAGE_mergedLAT-LATP_030816_rsID.bim"
#resultssubdir="${resultsdir}/results"
genopfx_altpop="${genepath}/${gene}_${altpop}"
rawpath_altpop="${genopfx_altpop}.raw"

# output files
#predictionfile="${outdir}/sage_predictions_${glmmethod}_${gene}.txt"
#lambdafile="${outdir}/sage_lambdas_${glmmethod}_${gene}.txt"
#weightsfile="${outdir}/sage_weights_${glmmethod}_${gene}.txt"
predictionfile="${resultssubdir}/geuvadis_predictions_${glmmethod}_${gene}.txt"
lambdafile="${resultssubdir}/geuvadis_lambdas_${glmmethod}_${gene}.txt"
weightsfile="${resultssubdir}/geuvadis_weights_${glmmethod}_${gene}.txt"
predictionfile_altpop="${resultssubdir}/geuvadis_predictinto_${altpop}_${glmmethod}_${gene}.txt"
predictionfile_samepop="${resultssubdir}/geuvadis_predictinto_${pop}_${glmmethod}_${gene}.txt"
h2filepfx_null="${resultssubdir}/geuvadis_h2_null_${pop}_${gene}"
h2filepfx="${resultssubdir}/geuvadis_h2_${pop}_${gene}"
hsqfile="${h2filepfx}.hsq"
hsqfile_null="${h2filepfx_null}.hsq"
h2file="${h2filepfx}.txt"
h2file_null="${h2filepfx_null}.txt"

echo "weightsfile = $weightsfile"
echo "predictionfile = $predictionfile"
echo "lambdafile = $lambdafile"
echo "predictionfile_altpop = $predictionfile_altpop"


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
 $PLINK \
    --bfile $bedfile \
    --chr ${chr} \
    --from-bp ${startpos} \
    --to-bp ${endpos} \
    --maf ${maf} \
    --hwe ${hwe} \
    --recode A \
    --make-bed \
    --out ${genopfx} \
    --threads ${nthreads} \
    --memory ${memory_limit_mb} \
    --keep ${subjectids} #\
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
 $Rscript $R_compute_new_weights ${rawpath} ${exprfile} ${gene} ${predictionfile} ${lambdafile} ${alpha} ${weightsfile} ${bimfile} ${nfolds}

# query return value of previous command
RETVAL=$?

# get list of SNPs to subset in alternate population
snps_to_extract="${genepath}/snps_to_extract_${altpop}_${gene}.txt"
cut -f 3 ${weightsfile} | grep -F "rs" | sort --temporary-directory $tmpdir | uniq > $snps_to_extract 


# create a PLINK RAW file
# this codes the dosage format required for glmnet
# here we use the genome-wide GEUVADIS genotype data with rsIDs
 $PLINK \
    --bfile $bedfile \
    --recode A \
    --make-bed \
    --out ${genopfx_altpop} \
    --threads ${nthreads} \
    --memory ${memory_limit_mb} \
    --keep ${subjectids_altpop} \
    --extract ${snps_to_extract} #\
#    --silent
#    --maf ${maf} \
#    --hwe ${hwe} \
#    --chr ${chr} \
#    --from-bp ${startpos} \
#    --to-bp ${endpos} \

# produce prediction from primary pop to alt pop
 $Rscript $R_predict_new_pop ${weightsfile} ${rawpath_altpop} ${predictionfile_altpop} ${gene} 

# query return value of previous command
let "RETVAL+=$?"

# also perform prediction from primary pop into itself
# different from true out-of-sample populations but systematically same as before
# want this to compare against quality of out-of-sample pops
 $Rscript $R_predict_new_pop ${weightsfile} ${rawpath} ${predictionfile_samepop} ${gene} 


# query return value of previous command
let "RETVAL+=$?"

# need cis heritability for this gene
# will run GCTA on the previous PLINK genotypes for current $pop 

# build GRM for the gene in $pop
$GCTA \
    --bfile ${genopfx} \
    --autosome \
    --make-grm-bin \
    --out ${genopfx} \
    --thread-num 1 # could make variable $nthreads = 1

# query return value of previous command
let "RETVAL+=$?"

# estimate null model 
$GCTA \
    --reml \
    --pheno ${phenofile} \
    --mpheno ${counter} \
    --out ${h2filepfx_null} \
    --thread-num 1 # $nthreads
    #--qcovar ${covar} \

# query return value of previous command
let "RETVAL+=$?"

# parse phenotyping variance for null model
echo -e "mygene\t$(grep "^Vp" $hsqfile_null | awk '{ print $2, "\t", $3 }')" > $h2file_null

# estimate genetic component of expression
reml_maxit=1000
reml_alg=2
maf_gcta="0.05"
$GCTA \
    --reml \
    --grm-bin ${genopfx} \
    --pheno ${phenofile} \
    --mpheno ${counter} \
    --out ${h2filepfx} \
    --reml-alg ${reml_alg} \
    --reml-maxit ${reml_maxit} \
    --maf ${maf_gcta} \
    --thread-num 1 # ${nthreads}
    #--grm-bin ${page_grm} \
    #--qcovar ${covar} \

# compute eigenvalues of GRM
# we need these to compute confidence intervals
# requires full spectral decomposition of GRM
# number of eigenvalues is lesser of (nsamples, nSNPs)
num_samples=$(wc -l ${genopfx}.fam | awk '{ print $1 }')
num_snps=$(wc -l ${genopfx}.bim | awk '{ print $1 }')
echo -e "Number of samples: ${num_samples}"
echo -e "Number of samples: ${num_snps}"
if [[ $num_samples -ge $num_snps ]]; then
    num_pcs=$num_snps
else
    num_pcs=$num_samples
fi
echo -e "Setting number of PCs from GCTA to $num_pcs"

$GCTA \
    --grm-bin ${genopfx} \
    --pca ${num_pcs} \
    --out ${genopfx} \
    --thread-num 1

# compute the CIs with ALBI/FIESTA
eigenvalues_file="${genopfx}.eigenval"
h2_CI_grid="${genopfx}.h2CI"
$PYTHON $FIESTA --kinship_eigenvalues ${eigenvalues_file} --estimate_grid 1000 --output_filename $h2_CI_grid

# parse results, inserting NA as necessary
# first parse null model, then parse genotype model 
pvar_null=$(grep "^Vp" $hsqfile_null | awk '{ print $2 }')
pse_null=$(grep "^Vp" $hsqfile_null | awk '{ print $3 }')
if [[ "${pvar_null}" = "" ]]; then pvar_null="NA"; fi
if [[ "${pse_null}" = "" ]]; then pse_null="NA"; fi
echo -e "$gene\t$pvar_null\t$pse_null" > $h2file_null

pvar=$(grep "^V(G)/Vp" $hsqfile | awk '{ print $2 }')
#pse=$(grep "^V(G)/Vp" $hsqfile | awk '{ print $3 }')
if [[ "${pvar}" = "" ]]; then pvar="NA"; fi
#if [[ "${pse}" = "" ]]; then pse="NA"; fi

# before writing nonnull results, parse the h2 CI results
# can do this with Rscript
# use h2 given from $pvar
# if it is NA then also call CIs as NA
# otherwise interpolate the CIs around h2
# do this by averaging the CIs corresponding to grid points just above and below h2
 #$Rscript -e "h2 = $pvar; if(is.na(h2)){ cat(\"Estimate\tCI_lower_bound\tCI_upper_bound\nNA\tNA\tNA\") } else { library(data.table); x = fread(\"$h2_CI_grid\"); apply(x[c(last(which(x\$Estimate < h2)), first(which(x\$Estimate > h2))),], 2, mean)}" > $h2_tempfile

h2_tempfile="${genepath}/${gene}.h2CIgrid"
if [[ "${pvar}" = "NA" ]]; then
    $Rscript -e "cat(\"Estimate\tCI_lower_bound\tCI_upper_bound\nNA\tNA\tNA\")" > $h2_tempfile
else
    $Rscript -e "library(data.table); h2 = $pvar; x = fread(\"$h2_CI_grid\"); apply(x[c(last(which(x\$Estimate < h2)), first(which(x\$Estimate > h2))),], 2, mean)" > $h2_tempfile
fi
h2_CI_lower=$(tail -n 1 $h2_tempfile | awk '{ print $2 }')
h2_CI_upper=$(tail -n 1 $h2_tempfile | awk '{ print $3 }')
#echo -e "$gene\t$pvar\t$pse" > $h2file
echo -e "$gene\t$pvar\t$h2_CI_lower\t$h2_CI_upper" > $h2file

# query return value of previous command
let "RETVAL+=$?"

#    # call R2 calculation script with default weights
#    #if [[ "$r2_default" != "0" -a "$r2_default" != ""]] ; then
#    if [[ ! -z "$r2_default" ]]; then
#        if [[ "$r2_default" != "0" ]]; then 
#            $Rscript $R_compute_r2 ${rawpath} ${exprfile} ${gene} ${ucsc_snpfile} ${r2resultsfile} "genotype_data"
#        fi
#    fi


# if return value is not 0, then previous command did not exit correctly
# create a status file notifying of error
# in contrary case, notify of success
if [ "$RETVAL" -ne "0" ];
then
  echo "ERROR" > ${logdir}/status.${glmmethod}.${gene}.${pop}
else
  echo "SUCCESS" > ${logdir}/status.${glmmethod}.${gene}.${pop}
fi

# DO NOT clean up scratch space, e.g., with `rm -rf ${scratchdir}`
# will subsequently compile results after this array job completes

#cat $predictionfile
#cat $lambdafile
#cat $weightsfile

# 
echo "Ended script at timestamp $(date)"

# ask job to report on itself
echo "job ${JOB_ID} report:"
qstat -j ${JOB_ID}
