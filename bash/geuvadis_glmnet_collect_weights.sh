#!/usr/bin/env bash                # -- what is the language of this shell
#$ -S /bin/bash                    # -- the shell for the job
##$ -M kevin.keys@ucsf.edu          # -- email status of this job to this address
##$ -m bes                          # -- email on beginning, end, and suspension of job
#$ -r y                            # -- tell the system that if a job crashes, it should be restarted
#$ -j y                            # -- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=1G                  # -- submits on nodes with enough free memory (required)
#$ -l arch=linux-x64               # -- SGE resources (CPU type)
#$ -l h_rt=0:10:00                 # -- runtime limit in hours 
#
# This script gathers results from a run of compute_new_predixcan_weights.sh.
# It produces two files:
# -- $resultsfile, which compiles prediction results for the current glmnet run;
# -- $weightsfile, which compiles nonzero prediction weights for the same glmnet run
#
# Call this script via qsub from compute_new_predixcan_weights.sh as follows:
#
# qsub -N glmnet.collect.results.${glmmethod} \
#      -o $logdir \
#      -e $logdir \
#      -v resultsfile=$resultsfile,weightfile=$weightfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir \
#      $glmnetdir/geuvadis_glmnet_collect_weights.sh
#
# coded by Kevin L. Keys (2018)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# user limits: -c max size of core files created
date
hostname
ulimit -c 0 

# script variables from qsub
weightsfile=$weightsfile
glmnetdir=$glmnetdir
glmmethod=$glmmethod
logdir=$logdir
lambdafile=$lambdafile
predictionfile=$predictionfile
outdir=$outdir
resultsdir=$resultsdir
resultssubdir=$resultssubdir
tmpdir=$tmpdir
pop=$pop
altpop=$altpop
predictionfile_altpop=$predictionfile_altpop
predictionfile_samepop=$predictionfile_samepop
predictionfile_header=$predictionfile_header

# check variables for debugging
#echo "resultsfile = $resultsfile"
#echo "weightsfile = $weightsfile"
#echo "glmnetdir = $glmnetdir"
#echo "glmmethod = $glmmethod"
#echo "logdir = $logdir"

## put header on prediction file
## then append results to file
#rm -f $predictionfile
#touch $predictionfile
#echo -e "$predictionfile_header" > $predictionfile
#cat ${resultssubdir}/geuvadis_predictions_${glmmethod}_ENSG* | sort --temporary-directory $tmpdir >> $predictionfile
#
## check previous command
#RETVAL=$?
#
## report on previous command
#echo "exit status after making prediction file: $RETVAL"
#
## clean up temporary directory
#rm -f $tmpdir/* 
#
## put header on lambda file
## then computed weights to file
## ensure when concatenating results that we discard the individual file headers
#rm -f $lambdafile
#touch $lambdafile
#echo -e "Gene\tHeld_out_Sample\tMean_MSE\tLambda\tPredicted_Expr" > $lambdafile
#cat ${resultssubdir}/geuvadis_lambdas_${glmmethod}_ENSG* | grep -v "Gene" | sort --temporary-directory $tmpdir >> $lambdafile
#
## check previous command
## compound with penultimate one
#let "RETVAL+=$?"
#echo "exit status after making prediction,lambda files: $RETVAL"
#
## clean up temporary directory
#rm -f $tmpdir/* 
#
## now compile weights file
#rm -f $weightsfile
#touch $weightsfile
#echo -e "Gene\tHeld_out_Sample\tSNP\tA1\tA2\tBeta" > $weightsfile
#cat ${resultssubdir}/geuvadis_weights_${glmmethod}_ENSG* | grep -v "Gene" | fgrep -v "NA"$'\t'"NA" | sort --temporary-directory $tmpdir >> $weightsfile
#
## report on previous commands
#let "RETVAL+=$?"
#echo "exit status after making prediction,lambda,weights files: $RETVAL"

# finally, compile external and internal prediction files
rm -f $predictionfile_altpop
touch $predictionfile_altpop
echo -e "Gene\tSubjectID\tPrediction_${pop}_into_${altpop}" > $predictionfile_altpop
cat ${resultssubdir}/geuvadis_predictinto_${altpop}_${glmmethod}_ENSG* | grep -v "SubjectID" | grep "ENS" | sort --temporary-directory $tmpdir >> $predictionfile_altpop

rm -f $predictionfile_samepop
touch $predictionfile_samepop
echo -e "Gene\tSubjectID\tPrediction_${pop}_into_${pop}" > $predictionfile_samepop
cat ${resultssubdir}/geuvadis_predictinto_${pop}_${glmmethod}_ENSG* | grep -v "SubjectID" | grep "ENS" | sort --temporary-directory $tmpdir >> $predictionfile_samepop


# report on previous commands
let "RETVAL+=$?"
echo "exit status after making prediction, lambda, weights, and external prediction files: $RETVAL"


# clean up temporary directory
rm -f $tmpdir/* 

# now it is safe to clean up scratch
rm -rf $outdir

if [ "$RETVAL" -ne "0" ];
then
    echo "ERROR" > ${logdir}/status.collectweights.${glmmethod}.${pop}
else
    echo "SUCCESS" > ${logdir}/status.collectweights.${glmmethod}.${pop}
fi
