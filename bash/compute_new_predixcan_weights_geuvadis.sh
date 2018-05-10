#!/usr/bin/env bash
#
# This script computes PrediXcan weights from gEUVADIS transcriptome data.
#
# Call:
#
# ./compute_new_predixcan_weights_geuvadis.sh $ALPHA
#
# where
# -- $ALPHA = "0.0", "0.5", or "1.0", used by glmnet to determine the training algorithm
# -- $DEFAULT_R2 is any nonempty nonzero value, used to compute R2/p-values with default GTEx weights (or not)
#
# coded by Kevin L. Keys (2017)

# terminate script on error
set -e

# parse command line arguments
alpha=$1
r2_default=$2

# set method based on alpha
# only checks the three admissible values for alpha
glmmethod=""
if [[ "$alpha" == "1.0" ]]; then
    glmmethod="lasso";
fi
if [[ "$alpha" == "0.5" ]]; then
    glmmethod="elasticnet";
fi
if [[ "$alpha" == "0.0" ]]; then
    glmmethod="ridge";
fi

# make sure that glmmethod was set
# $glmmethod is empty if alpha was not given an acceptable value
if [[ "$glmmethod" == "" ]]; then
    echo -e "usage:\n\tsource compute_new_predixcan_weights.sh \$ALPHA [\$DEFAULT]\n"
    echo -e "where \$ALPHA = {0.0, 0.5, 1.0} (required)\nand DEFAULT != 0 (optional, if calulation with default weights desired)\n"
    return 1;
fi

# script static directories
MYHOME="/netapp/home/kkeys"
rnaseqdir="${MYHOME}/gala_sage/rnaseq"
geuvadisdir="${rnaseqdir}/geuvadis"
gctadir="${rnaseqdir}/gcta"
bindir="${MYHOME}/bin"
datadir="${rnaseqdir}/data"
glmnetdir="${rnaseqdir}/glmnet"
genotypedir="${MYHOME}/gala_sage/genotypes/mergedLAT-LATP/SAGE_merge"
logdir="${glmnetdir}/log"
#outdir="${glmnetdir}/${glmmethod}/genes"
outdir_eur="/scratch/kkeys/${glmmethod}/genes/eur373"
outdir_yri="/scratch/kkeys/${glmmethod}/genes/yri89"
#imputegenodir="${MYHOME}/gala_sage/genotypes/SAGE/IMPUTE_HRC_1.1_PLINK"
imputegenodir="${MYHOME}/gala_sage/genotypes/gEUVADIS"
codedir="${MYHOME}/gala_sage/code"
resultsdir_eur="${glmnetdir}/${glmmethod}/eur373"
resultsdir_yri="${glmnetdir}/${glmmethod}/yri89"
resultssubdir_eur="${resultsdir_eur}/results"
resultssubdir_yri="${resultsdir_yri}/results"
tmpdir="${MYHOME}/tmp"

# make output and results directories, in case they doesn't exist
mkdir -p $outdir_eur
mkdir -p $outdir_yri
mkdir -p $resultsdir_eur
mkdir -p $resultsdir_eur
mkdir -p $resultssubdir_eur
mkdir -p $resultssubdir_yri
mkdir -p $tmpdir

# script files
exprfile_eur="${geuvadisdir}/geuvadis.eur373.RPKM.invnorm.txt"
exprfile_yri="${geuvadisdir}/geuvadis.yri89.RPKM.invnorm.txt"
#predictionfile="${glmnetdir}/sage_${glmmethod}_predictions.txt"
#lambdafile="${glmnetdir}/sage_${glmmethod}_lambdas.txt"
#weightsfile="${glmnetdir}/sage_${glmmethod}_weights.txt"
ucsc_snpfile="${glmnetdir}/sage_rnaseq_UCSC_allsnps_merged_wgs_rsid.txt"
r2resultsfile="${glmnetdir}/sage_gtex_r2_defaultweights.txt"
#sagegenopfx="${genotypedir}/SAGE_mergedLAT-LATP_030816_rsID"
#genelist="${gctadir}/genelist_plusmin1Mb.txt"
genelist="${geuvadisdir}/human_ens_GRCh37_genechrpos_plusmin500kb.txt"
#genelist="${geuvadisdir}/human_ens_GRCh37_genechrpos_plusmin500kb_failedgenes.txt"
#sage_lab2nwd="${datadir}/sage_lab2nwd.txt"
#subjectids="${datadir}/sage_39_subjectids.txt"
#nwdids="${datadir}/sage_39_nwdids.txt"
subjectids_eur="${geuvadisdir}/geuvadis.eur373.sampleids.txt"
subjectids_yri="${geuvadisdir}/geuvadis.yri89.sampleids.txt"
predictionfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictions.txt"
predictionfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictions.txt"
predictionfile_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89.txt"
predictionfile_eur2eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_eur373.txt"
predictionfile_yri2eur="${resultsdir_eur}/geuvadis_${glmmethod}_yri89_predictinto_eur373.txt"
predictionfile_yri2yri="${resultsdir_eur}/geuvadis_${glmmethod}_yri89_predictinto_yri89.txt"
lambdafile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lambdas.txt"
lambdafile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lambdas.txt"
weightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights.txt"
weightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights.txt"
num_pred_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_numpred.txt"
num_pred_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_numpred.txt"
newweightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights_noNA.txt"
newweightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights_noNA.txt"
out_lm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lm_predvmeas_results.txt"
out_lm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lm_predvmeas_results.txt"
out_genelm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_genelm_predvmeas_results.txt"
out_genelm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_genelm_predvmeas_results.txt"
out_lm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_lm_predvmeas_results.txt"
out_lm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_lm_predvmeas_results.txt"
out_genelm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_genelm_predvmeas_results.txt"
out_genelm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_genelm_predvmeas_results.txt"

# script locations
R_compute_new_weights="${codedir}/geuvadis_gtex_compute_new_weights.R"
R_compute_r2="${codedir}/geuvadis_gtex_compute_r2.R"
R_glmnet_postprocess="${codedir}/geuvadis_glmnet_postprocess.R"
R_predict_new_pop="${codedir}/geuvadis_predict_in_altpop.R"
BASH_run_glmnet="${codedir}/qsub_run_geuvadis_glmnet.sh"
BASH_collect_weights="${codedir}/geuvadis_glmnet_collect_weights.sh"
BASH_postprocess_glmnet="${codedir}/qsub_sage_glmnet_postprocess.sh"

# script binaries
PLINK="${bindir}/plink"
Rscript="${bindir}/Rscript"
#Rscript="/usr/bin/Rscript"

# script variables
maf="0.01"
hwe="0.0001"
geno="0.03"
nthreads=1
memory_limit="2G"
memory_limit_mb="2000" # manually coordinate this with $memory_limit!!!
scratch_memory="2G"
#h_rt="12:00:00"
h_rt="06:00:00"
discard_ratio="0.5" # desired min ratio of samples with nonmissing LOOCV predictions, used in postprocessing, e.g. 0.5 = 50%
nsamples_eur="373"
nsamples_yri="89"

# need headers for prediction files
predictionfile_header_eur="Gene\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101\tHG00102\tHG00103\tHG00104\tHG00105\tHG00106\tHG00108\tHG00109\tHG00110\tHG00111\tHG00112\tHG00114\tHG00115\tHG00116\tHG00117\tHG00118\tHG00119\tHG00120\tHG00121\tHG00122\tHG00123\tHG00124\tHG00125\tHG00126\tHG00127\tHG00128\tHG00129\tHG00130\tHG00131\tHG00132\tHG00133\tHG00134\tHG00135\tHG00136\tHG00137\tHG00138\tHG00139\tHG00141\tHG00142\tHG00143\tHG00145\tHG00146\tHG00148\tHG00149\tHG00150\tHG00151\tHG00152\tHG00154\tHG00155\tHG00156\tHG00157\tHG00158\tHG00159\tHG00160\tHG00171\tHG00173\tHG00174\tHG00176\tHG00177\tHG00178\tHG00179\tHG00180\tHG00181\tHG00182\tHG00183\tHG00185\tHG00186\tHG00187\tHG00188\tHG00189\tHG00231\tHG00232\tHG00233\tHG00234\tHG00235\tHG00236\tHG00238\tHG00239\tHG00240\tHG00242\tHG00243\tHG00244\tHG00245\tHG00246\tHG00247\tHG00249\tHG00250\tHG00251\tHG00252\tHG00253\tHG00255\tHG00256\tHG00257\tHG00258\tHG00259\tHG00260\tHG00261\tHG00262\tHG00263\tHG00264\tHG00265\tHG00266\tHG00267\tHG00268\tHG00269\tHG00271\tHG00272\tHG00273\tHG00274\tHG00275\tHG00276\tHG00277\tHG00278\tHG00280\tHG00281\tHG00282\tHG00284\tHG00285\tHG00306\tHG00308\tHG00309\tHG00310\tHG00311\tHG00312\tHG00313\tHG00315\tHG00319\tHG00320\tHG00321\tHG00323\tHG00324\tHG00325\tHG00326\tHG00327\tHG00328\tHG00329\tHG00330\tHG00331\tHG00332\tHG00334\tHG00335\tHG00336\tHG00337\tHG00338\tHG00339\tHG00341\tHG00342\tHG00343\tHG00344\tHG00345\tHG00346\tHG00349\tHG00350\tHG00351\tHG00353\tHG00355\tHG00356\tHG00358\tHG00359\tHG00360\tHG00361\tHG00362\tHG00364\tHG00365\tHG00366\tHG00367\tHG00369\tHG00371\tHG00372\tHG00373\tHG00375\tHG00376\tHG00377\tHG00378\tHG00379\tHG00380\tHG00381\tHG00382\tHG00383\tHG00384\tHG01334\tHG01789\tHG01790\tHG01791\tHG02215\tNA06984\tNA06985\tNA06986\tNA06989\tNA06994\tNA07037\tNA07048\tNA07051\tNA07056\tNA07346\tNA07347\tNA07357\tNA10847\tNA10851\tNA11829\tNA11830\tNA11831\tNA11832\tNA11840\tNA11843\tNA11881\tNA11892\tNA11893\tNA11894\tNA11918\tNA11920\tNA11930\tNA11931\tNA11992\tNA11993\tNA11994\tNA11995\tNA12004\tNA12005\tNA12006\tNA12043\tNA12044\tNA12045\tNA12058\tNA12144\tNA12154\tNA12155\tNA12156\tNA12234\tNA12249\tNA12272\tNA12273\tNA12275\tNA12282\tNA12283\tNA12286\tNA12287\tNA12340\tNA12341\tNA12342\tNA12347\tNA12348\tNA12383\tNA12399\tNA12400\tNA12413\tNA12489\tNA12546\tNA12716\tNA12717\tNA12718\tNA12749\tNA12750\tNA12751\tNA12760\tNA12761\tNA12762\tNA12763\tNA12775\tNA12776\tNA12777\tNA12778\tNA12812\tNA12813\tNA12814\tNA12815\tNA12827\tNA12829\tNA12830\tNA12842\tNA12843\tNA12872\tNA12873\tNA12874\tNA12889\tNA12890\tNA20502\tNA20503\tNA20504\tNA20505\tNA20506\tNA20507\tNA20508\tNA20509\tNA20510\tNA20512\tNA20513\tNA20514\tNA20515\tNA20516\tNA20517\tNA20518\tNA20519\tNA20520\tNA20521\tNA20524\tNA20525\tNA20527\tNA20528\tNA20529\tNA20530\tNA20531\tNA20532\tNA20534\tNA20535\tNA20536\tNA20537\tNA20538\tNA20539\tNA20540\tNA20541\tNA20542\tNA20543\tNA20544\tNA20581\tNA20582\tNA20585\tNA20586\tNA20588\tNA20589\tNA20752\tNA20754\tNA20756\tNA20757\tNA20758\tNA20759\tNA20760\tNA20761\tNA20765\tNA20766\tNA20768\tNA20769\tNA20770\tNA20771\tNA20772\tNA20773\tNA20774\tNA20778\tNA20783\tNA20785\tNA20786\tNA20787\tNA20790\tNA20792\tNA20795\tNA20796\tNA20797\tNA20798\tNA20799\tNA20800\tNA20801\tNA20802\tNA20803\tNA20804\tNA20805\tNA20806\tNA20807\tNA20808\tNA20809\tNA20810\tNA20811\tNA20812\tNA20813\tNA20814\tNA20815\tNA20816\tNA20819\tNA20826\tNA20828" 
predictionfile_header_yri=""

### switches to control which parts of script to run
### set to 0 in order to run
run_compute_weights_eur=0
run_compute_weights_yri=1

# -------------------- #
# start script
echo "Start time: $(date)"

# -------------------- #
# make necessary output directories
if [[ ! -d "$logdir" ]]; then mkdir -p $logdir; fi
#if [[ ! -d "$outdir" ]]; then mkdir -p $outdir; fi


# -------------------- #
## do the same for default weights file
#if [[ ! -z "$r2_default" ]]; then
#    if [[ "$r2_default" != "0" ]]; then 
#        touch $r2resultsfile
#        echo -e "gene\tn.snps\tR2\tpval" > $r2resultsfile
#    fi
#fi



# -------------------- #
# compute new weights for each gene
# this entails running glmnet on a PLINK RAW genotype dosage file
# and extracting the expression level from the UCSC BED file stored at $exprfile

## how many genes are in this list?
## remember that this file has a header so line count is off by 1 
nGenes=$(wc -l $genelist | cut -f 1 -d " ")
###let "nGenes=nGenes-1"

### for testing
#nGenes=1


## -N glmnet.${gene}    ##--- wait until array job lowComplex_byPop.${PR,AF,MX}
## -e $glmnetdir/log    ##--- where to put error files
## -o $glmnetdir/log    ##--- where to put output files
## -v ...               ##--- variables to pass to script
## -t ...               ##--- array job components, one for each gene
## -l mem_free=1G       ## -- submits on nodes with enough free memory, e.g. 1 gigabyte (required)
## -l h_rt=0:29:59      ## -- runtime limit in hours 
 
### first do EUR
if [[ "${run_compute_weights_eur}" -eq "0" ]]; then

    # set variables to population-specific parameters
    # point prediction, weight, lambda, expression, subject files to EUR
    # set alternative population parameters to YRI

    h_rt="00:29:59" ### REVISE ME WHEN RETRAINING WEIGHTS!!!!
    pop="eur373" 
    altpop="yri89"
    subjectids=$subjectids_eur
    exprfile=$exprfile_eur
    outdir=$outdir_eur
    resultssubdir=$resultssubdir_eur
    resultsdir=$resultsdir_eur
    subjectids_altpop=$subjectids_yri
    predictionfile=$predictionfile_eur
    lambdafile=$lambdafile_eur
    predictionfile_altpop=$predictionfile_eur2yri
    predictionfile_samepop=$predictionfile_eur2eur
    altpop_out_lm_file=$out_lm_file_eur2yri
    altpop_out_genelm_file=$out_genelm_file_eur2yri
    altpop_exprfile=$exprfile_yri
    predictionfile_header=$predictionfile_header_eur
    weightsfile=$weightsfile_eur
    newweightsfile=$newweightsfile_eur
    out_genelm_file=$out_genelm_file_eur
    num_pred_file=$num_pred_file_eur
    nsamples=$nsamples_eur
    out_lm_file=$out_lm_file_eur

    qsub -N glmnet.${glmmethod}.${pop} \
         -v genelist=$genelist,subjectids=$subjectids,Rscript=$Rscript,R_compute_new_weights=$R_compute_new_weights,exprfile=$exprfile,logdir=$logdir,alpha=$alpha,gctadir=$gctadir,glmmethod=$glmmethod,outdir=$outdir,imputegenodir=$imputegenodir,maf=$maf,hwe=$hwe,nthreads=$nthreads,PLINK=$PLINK,memory_limit_mb=$memory_limit_mb,resultssubdir=$resultssubdir,resultsdir=$resultsdir,pop=$pop,altpop=$altpop,tmpdir=$tmpdir,subjectids_altpop=$subjectids_altpop,R_predict_new_pop=$R_predict_new_pop \
         -t 1-$nGenes \
         -e $logdir \
         -o $logdir \
         -l mem_free=$memory_limit \
         -l scratch=$scratch_memory \
         -l h_rt=$h_rt \
         $BASH_run_glmnet

    # beefy calculations done, don't need as much compute time
    # set $h_rt to < 30 min to schedule in QB3 speedy queue 
    h_rt="00:29:59"
    h_rt="06:00:00"

    # collect output
    qsub -N glmnet.${glmmethod}.collect.weights.${pop} \
         -hold_jid glmnet.${glmmethod}.${pop} \
         -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile,lambdafile=$lambdafile,outdir=$outdir,resultsdir=$resultsdir,resultssubdir=$resultssubdir,tmpdir=$tmpdir,pop=$pop,altpop=$altpop,predictionfile_altpop=$predictionfile_altpop,predictionfile_samepop=$predictionfile_samepop,predictionfile_header=$predictionfile_header \
         -o $logdir \
         -e $logdir \
         -l mem_free=$memory_limit \
         -l h_rt=$h_rt \
         $BASH_collect_weights

    # process output file
    qsub -N glmnet.${glmmethod}.postprocess.${pop} \
        -hold_jid glmnet.${glmmethod}.${pop},glmnet.${glmmethod}.collect.weights.${pop} \
        -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,newweightsfile=$newweightsfile,Rscript=$Rscript,R_glmnet_postprocess=$R_glmnet_postprocess,discard_ratio=$discard_ratio,num_pred_file=$num_pred_file,nsamples=$nsamples,predictionfile=$predictionfile,exprfile=$exprfile,out_lm_file=$out_lm_file,out_genelm_file=$out_genelm_file,logdir=$logdir,pop=$pop,predictionfile_altpop=$predictionfile_altpop,altpop_exprfile=$altpop_exprfile,altpop_out_lm_file=$altpop_out_lm_file,altpop_out_genelm_file=$altpop_out_genelm_file \
        -o $logdir \
        -e $logdir \
         -l mem_free=$memory_limit \
         -l h_rt=$h_rt \
        ${codedir}/qsub_sage_glmnet_postprocess.sh
fi

### now do YRI
################   2018-05-07: NEED TO UPDATE!!!!!
if [[ "${run_compute_weights_yri}" -eq "0" ]]; then

    pop="yri89"
    altpop="eur373"
    subjectids_altpop=$subjectids_eur

    # smaller sample size --> less compute time necessary
    h_rt="00:29:59"
    qsub -N glmnet.${glmmethod}.${pop} \
         -v genelist=$genelist,subjectids=$subjectids_yri,Rscript=$Rscript,R_compute_new_weights=$R_compute_new_weights,exprfile=$exprfile_yri,logdir=$logdir,alpha=$alpha,gctadir=$gctadir,glmmethod=$glmmethod,outdir=$outdir_yri,imputegenodir=$imputegenodir,maf=$maf,hwe=$hwe,nthreads=$nthreads,PLINK=$PLINK,memory_limit_mb=$memory_limit_mb,resultssubdir=$resultssubdir_yri,resultsdir=$resultsdir_yri,pop=$pop \
         -t 1-$nGenes \
         -e $logdir \
         -o $logdir \
         -l mem_free=$memory_limit \
         -l scratch=$scratch_memory \
         -l h_rt=$h_rt \
         ${codedir}/qsub_run_geuvadis_glmnet.sh

    qsub -N glmnet.${glmmethod}.collect.weights.${pop} \
         -hold_jid glmnet.${glmmethod}.${pop} \
         -v weightsfile=$weightsfile_yri,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile_yri,lambdafile=$lambdafile_yri,outdir=$outdir_yri,resultsdir=$resultsdir_yri,resultssubdir=$resultssubdir_yri \
         -o $logdir \
         -e $logdir \
         -l mem_free=$memory_limit \
         -l h_rt=$h_rt \
         ${codedir}/geuvadis_glmnet_collect_weights.sh


    qsub -N glmnet.${glmmethod}.postprocess.${pop} \
        -hold_jid glmnet.${glmmethod}.${pop},glmnet.${glmmethod}.collect.weights.${pop} \
        -v weightsfile=$weightsfile_yri,glmnetdir=$glmnetdir,newweightsfile=$newweightsfile_yri,Rscript=$Rscript,R_glmnet_postprocess=$R_glmnet_postprocess,discard_ratio=$discard_ratio,num_pred_file=$num_pred_file_yri,nsamples=$nsamples_yri,predictionfile=$predictionfile_yri,exprfile=$exprfile_yri,out_lm_file=$out_lm_file_yri,out_genelm_file=$out_genelm_file_yri,logdir=$logdir \
        -o $logdir \
        -e $logdir \
         -l mem_free=$memory_limit \
         -l h_rt=$h_rt \
        ${codedir}/qsub_geuvadis_glmnet_postprocess.sh
fi

# end script
echo "End time: $(date)"
