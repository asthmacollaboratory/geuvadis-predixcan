#!/usr/bin/env bash
# ==============================================================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script is called from the main GEUVADIS prediction model training pipeline. IT DOES NOT WORK OUTSIDE OF THE SCRIPT.
# The requisite variables in the following qsub arguments are set in the pipeline.
# This script schedules three jobs:
#     1. compute new weights from a training set of genotypes and gene expression measures
#     2. collect the individual gene files and compile them into a handful of files
#     3. postprocess the results, generate imputation R2 and correlations
# ==============================================================================================================================

# estimate new prediction weights with glmnet, parallelized across genes
qsub -N glmnet.${glmmethod}.${pop} \
     -v genelist=$genelist,subjectids=$subjectids,Rscript=$Rscript,R_compute_new_weights=$R_compute_new_weights,exprfile=$exprfile,logdir=$logdir,alpha=$alpha,gctadir=$gctadir,glmmethod=$glmmethod,outdir=$outdir,imputegenodir=$imputegenodir,maf=$maf,hwe=$hwe,nthreads=$nthreads,PLINK=$PLINK,memory_limit_mb=$memory_limit_mb,resultssubdir=$resultssubdir,resultsdir=$resultsdir,pop=$pop,altpop=$altpop,tmpdir=$tmpdir,subjectids_altpop=$subjectids_altpop,R_predict_new_pop=$R_predict_new_pop,GCTA=$GCTA,phenofile=$phenofile,FIESTA=$PYTHON_fiesta,PYTHON=$PYTHON,nfolds=$nfolds \
     -t 1-$nGenes \
     -e $logdir \
     -o $logdir \
     -l mem_free=$memory_limit \
     -l scratch=$scratch_memory \
     -l h_rt=$h_rt \
     $BASH_run_glmnet

# beefy calculations done, may not need as much compute time
# can set $h_rt to < 30 min to schedule in QB3 speedy queue 
# nota bene: jobs in short.q killed at 30 min regardless of completion status!!!
#h_rt="00:29:59"
h_rt="06:00:00"

# -------------------- #
# for each gene, collect weights, lambdas, predictions, etc.
# merge these results into single file for weights, single file for lambdas, etc.
qsub -N glmnet.${glmmethod}.collect.weights.${pop} \
     -hold_jid glmnet.${glmmethod}.${pop} \
     -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile,lambdafile=$lambdafile,outdir=$outdir,resultsdir=$resultsdir,resultssubdir=$resultssubdir,tmpdir=$tmpdir,pop=$pop,altpop=$altpop,predictionfile_altpop=$predictionfile_altpop,predictionfile_samepop=$predictionfile_samepop,predictionfile_header="${predictionfile_header}",h2file=$h2file,h2file_null=$h2file_null \
     -o $logdir \
     -e $logdir \
     -l mem_free=$memory_limit \
     -l h_rt=$h_rt \
     $BASH_collect_weights

# -------------------- #
# process output files
qsub -N glmnet.${glmmethod}.postprocess.${pop} \
    -hold_jid glmnet.${glmmethod}.${pop},glmnet.${glmmethod}.collect.weights.${pop} \
    -v glmmethod=$glmmethod,weightsfile=$weightsfile,glmnetdir=$glmnetdir,newweightsfile=$newweightsfile,Rscript=$Rscript,R_glmnet_postprocess=$R_glmnet_postprocess,discard_ratio=$discard_ratio,num_pred_file=$num_pred_file,nsamples=$nsamples,predictionfile=$predictionfile,exprfile=$exprfile,out_lm_file=$out_lm_file,out_genelm_file=$out_genelm_file,logdir=$logdir,pop=$pop,predictionfile_altpop=$predictionfile_altpop,altpop_exprfile=$altpop_exprfile,altpop_out_lm_file=$altpop_out_lm_file,altpop_out_genelm_file=$altpop_out_genelm_file \
    -o $logdir \
    -e $logdir \
     -l mem_free=$memory_limit \
     -l h_rt=$h_rt \
     $BASH_postprocess_glmnet
