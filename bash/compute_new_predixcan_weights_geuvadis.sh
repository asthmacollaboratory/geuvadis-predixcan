#!/usr/bin/env bash
# ==============================================================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script computes PrediXcan weights from gEUVADIS transcriptome data.
#
# Call:
#
# ./compute_new_predixcan_weights_geuvadis.sh $ALPHA
#
# where
# -- $ALPHA = "0.0", "0.5", or "1.0", used by glmnet to determine the training algorithm
# ==============================================================================================================================

# ==============================================================================================================================
# parse command line arguments
# ==============================================================================================================================

# terminate script on error
set -e

# parse command line arguments
alpha=$1

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
    echo -e "usage:\n\tsource compute_new_predixcan_weights.sh \$ALPHA\n"
    echo -e "where \$ALPHA = {0.0, 0.5, 1.0} (required)\n"
    return 1;
fi

# ==============================================================================================================================
# set script variables
# ==============================================================================================================================

# script static directories
MYHOME="/netapp/home/kkeys"
rnaseqdir="${MYHOME}/gala_sage/rnaseq"
geuvadisdir="${rnaseqdir}/geuvadis"
gctadir="${rnaseqdir}/gcta"
bindir="${MYHOME}/bin"
datadir="${rnaseqdir}/data"
glmnetdir="${rnaseqdir}/glmnet"
logdir="${glmnetdir}/log"
codedir="${MYHOME}/gala_sage/code"
imputegenodir="${MYHOME}/gala_sage/genotypes/gEUVADIS"
outdir_eur="/scratch/kkeys/${glmmethod}/genes/eur373"
outdir_eur278="/scratch/kkeys/${glmmethod}/genes/eur278"
outdir_yri="/scratch/kkeys/${glmmethod}/genes/yri89"
outdir_tsi="/scratch/kkeys/${glmmethod}/genes/tsi"
outdir_gbr="/scratch/kkeys/${glmmethod}/genes/gbr"
outdir_fin="/scratch/kkeys/${glmmethod}/genes/fin"
outdir_ceu="/scratch/kkeys/${glmmethod}/genes/ceu"
resultsdir_eur="${glmnetdir}/${glmmethod}/eur373"
resultsdir_eur278="${glmnetdir}/${glmmethod}/eur278"
resultsdir_yri="${glmnetdir}/${glmmethod}/yri89"
resultsdir_ceu="${glmnetdir}/${glmmethod}/ceu92"
resultsdir_gbr="${glmnetdir}/${glmmethod}/gbr96"
resultsdir_fin="${glmnetdir}/${glmmethod}/fin95"
resultsdir_tsi="${glmnetdir}/${glmmethod}/tsi93"
resultsdir_crosspop="${glmnetdir}/${glmmethod}/crosspop"


resultsdir_ceu89="${resultsdir_crosspop}/ceu89"
resultsdir_tsi89="${resultsdir_crosspop}/tsi89"
resultsdir_gbr89="${resultsdir_crosspop}/gbr89"
resultsdir_fin89="${resultsdir_crosspop}/fin89"
resultssubdir_eur="${resultsdir_eur}/results"
resultssubdir_eur278="${resultsdir_eur278}/results"
resultssubdir_yri="${resultsdir_yri}/results"
resultssubdir_ceu="${resultsdir_ceu}/results"
resultssubdir_gbr="${resultsdir_gbr}/results"
resultssubdir_fin="${resultsdir_fin}/results"
resultssubdir_tsi="${resultsdir_tsi}/results"
resultssubdir_ceu89="${resultsdir_ceu89}/results"
resultssubdir_gbr89="${resultsdir_gbr89}/results"
resultssubdir_fin89="${resultsdir_fin89}/results"
resultssubdir_tsii89="${resultsdir_tsi89}/results"
tmpdir="${MYHOME}/tmp"
crosspop_dir="${geuvadisdir}/crosspop"

# make output and results directories, in case they doesn't exist
mkdir -p $outdir_eur
mkdir -p $outdir_eur278
mkdir -p $outdir_yri
mkdir -p $resultsdir_eur
mkdir -p $resultsdir_eur278
mkdir -p $resultsdir_yri
mkdir -p $resultsdir_ceu
mkdir -p $resultsdir_gbr
mkdir -p $resultsdir_fin
mkdir -p $resultsdir_tsi
mkdir -p $resultssubdir_eur
mkdir -p $resultssubdir_eur278
mkdir -p $resultssubdir_yri
mkdir -p $resultssubdir_ceu
mkdir -p $resultssubdir_gbr
mkdir -p $resultssubdir_fin
mkdir -p $resultssubdir_tsi
mkdir -p $tmpdir
mkdir -p $crosspop_dir

# script files
exprfile_eur="${geuvadisdir}/geuvadis.eur373.RPKM.invnorm.txt"
exprfile_eur278="${geuvadisdir}/geuvadis.eur278.RPKM.invnorm.txt"
exprfile_yri="${geuvadisdir}/geuvadis.yri89.RPKM.invnorm.txt"
exprfile_fin="${geuvadisdir}/geuvadis.fin95.RPKM.invnorm.txt"
phenofile_eur="${geuvadisdir}/geuvadis.eur373.RPKM.invnorm.pheno"
phenofile_eur278="${geuvadisdir}/geuvadis.eur278.RPKM.invnorm.pheno"
phenofile_yri="${geuvadisdir}/geuvadis.yri89.RPKM.invnorm.pheno"
phenofile_fin="${geuvadisdir}/geuvadis.fin95.RPKM.invnorm.pheno"
genelist="${geuvadisdir}/human_ens_GRCh37_genechrpos_plusmin500kb.txt"
subjectids_eur="${geuvadisdir}/geuvadis.eur373.sampleids.txt"
subjectids_yri="${geuvadisdir}/geuvadis.yri89.sampleids.txt"
subjectids_fin="${geuvadisdir}/geuvadis.fin95.sampleids.txt"
subjectids_eur278="${geuvadisdir}/geuvadis.eur278.sampleids.txt"
predictionfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictions.txt"
predictionfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictions.txt"
predictionfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictions.txt"
predictionfile_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89.txt"
predictionfile_eur2eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_eur373.txt"
predictionfile_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89.txt"
predictionfile_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95.txt"
predictionfile_eur278toeur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_eur278.txt"
predictionfile_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373.txt"
predictionfile_yri2yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_yri89.txt"
lambdafile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lambdas.txt"
lambdafile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_lambdas.txt"
lambdafile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lambdas.txt"
weightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights.txt"
weightsfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_weights.txt"
weightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights.txt"
num_pred_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_numpred.txt"
num_pred_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_numpred.txt"
num_pred_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_numpred.txt"
newweightsfile_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_weights_noNA.txt"
newweightsfile_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_weights_noNA.txt"
newweightsfile_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_weights_noNA.txt"
out_lm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_lm_predvmeas_results.txt"
out_lm_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_lm_predvmeas_results.txt"
out_lm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_lm_predvmeas_results.txt"
out_genelm_file_eur="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_genelm_predvmeas_results.txt"
out_genelm_file_eur278="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_genelm_predvmeas_results.txt"
out_genelm_file_yri="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_genelm_predvmeas_results.txt"
out_lm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_lm_predvmeas_results.txt"
out_lm_file_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95_lm_predvmeas_results.txt"
out_lm_file_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89_lm_predvmeas_results.txt"
out_lm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_lm_predvmeas_results.txt"
out_genelm_file_eur2yri="${resultsdir_eur}/geuvadis_${glmmethod}_eur373_predictinto_yri89_genelm_predvmeas_results.txt"
out_genelm_file_eur278toyri="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_yri89_genelm_predvmeas_results.txt"
out_genelm_file_eur278tofin="${resultsdir_eur278}/geuvadis_${glmmethod}_eur278_predictinto_fin95_genelm_predvmeas_results.txt"
out_genelm_file_yri2eur="${resultsdir_yri}/geuvadis_${glmmethod}_yri89_predictinto_eur373_genelm_predvmeas_results.txt"
h2file_eur="${resultsdir_eur}/geuvadis_h2_eur373.txt"
h2file_null_eur="${resultsdir_eur}/geuvadis_h2_null_eur373.txt"
h2file_eur278="${resultsdir_eur278}/geuvadis_h2_eur278.txt"
h2file_null_eur278="${resultsdir_eur278}/geuvadis_h2_null_eur278.txt"
h2file_yri="${resultsdir_yri}/geuvadis_h2_yri89.txt"
h2file_null_yri="${resultsdir_yri}/geuvadis_h2_null_yri89.txt"

# script locations
R_compute_new_weights="${codedir}/geuvadis_glmnet_compute_new_weights.R"
R_compute_r2="${codedir}/geuvadis_gtex_compute_r2.R"
R_glmnet_postprocess="${codedir}/geuvadis_glmnet_postprocess.R"
R_subsample_eur="${codedir}/geuvadis_subsample_eur373.R" 
R_predict_new_pop="${codedir}/geuvadis_predict_in_altpop.R"
R_subsample_pop="${codedir}/geuvadis_subsample_onepop.R" 
BASH_run_glmnet="${codedir}/qsub_geuvadis_run_glmnet.sh"
BASH_collect_weights="${codedir}/qsub_geuvadis_glmnet_collect_weights.sh"
BASH_postprocess_glmnet="${codedir}/qsub_geuvadis_glmnet_postprocess.sh"
PYTHON_fiesta="${HOME}/git/albi/fiesta.py"

# script binaries
PLINK="${bindir}/plink"
Rscript="${bindir}/Rscript"
GCTA="${bindir}/gcta64"
PYTHON="${bindir}/python2"

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
nsamples_eur278="278"
nsamples_yri="89"
nfolds_eur="10"
nfolds_eur278="10"
nfolds_yri="89"

### parameters for EUR89
# set variables to population-specific parameters
# point prediction, weight, lambda, expression, subject files to EUR
# set alternative population parameters to YRI89 
h_rt="23:59:59"
pop="eur89" 
nresample="100"
eur89_dir="${geuvadisdir}/eur89"
nfolds_eur89=${nfolds_yri} ## ensure that the fold number matches the number used for training AFR 
nsamples_eur89="89"

# need headers for prediction files
predictionfile_header_eur373="Gene\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101\tHG00102\tHG00103\tHG00104\tHG00105\tHG00106\tHG00108\tHG00109\tHG00110\tHG00111\tHG00112\tHG00114\tHG00115\tHG00116\tHG00117\tHG00118\tHG00119\tHG00120\tHG00121\tHG00122\tHG00123\tHG00124\tHG00125\tHG00126\tHG00127\tHG00128\tHG00129\tHG00130\tHG00131\tHG00132\tHG00133\tHG00134\tHG00135\tHG00136\tHG00137\tHG00138\tHG00139\tHG00141\tHG00142\tHG00143\tHG00145\tHG00146\tHG00148\tHG00149\tHG00150\tHG00151\tHG00152\tHG00154\tHG00155\tHG00156\tHG00157\tHG00158\tHG00159\tHG00160\tHG00171\tHG00173\tHG00174\tHG00176\tHG00177\tHG00178\tHG00179\tHG00180\tHG00181\tHG00182\tHG00183\tHG00185\tHG00186\tHG00187\tHG00188\tHG00189\tHG00231\tHG00232\tHG00233\tHG00234\tHG00235\tHG00236\tHG00238\tHG00239\tHG00240\tHG00242\tHG00243\tHG00244\tHG00245\tHG00246\tHG00247\tHG00249\tHG00250\tHG00251\tHG00252\tHG00253\tHG00255\tHG00256\tHG00257\tHG00258\tHG00259\tHG00260\tHG00261\tHG00262\tHG00263\tHG00264\tHG00265\tHG00266\tHG00267\tHG00268\tHG00269\tHG00271\tHG00272\tHG00273\tHG00274\tHG00275\tHG00276\tHG00277\tHG00278\tHG00280\tHG00281\tHG00282\tHG00284\tHG00285\tHG00306\tHG00308\tHG00309\tHG00310\tHG00311\tHG00312\tHG00313\tHG00315\tHG00319\tHG00320\tHG00321\tHG00323\tHG00324\tHG00325\tHG00326\tHG00327\tHG00328\tHG00329\tHG00330\tHG00331\tHG00332\tHG00334\tHG00335\tHG00336\tHG00337\tHG00338\tHG00339\tHG00341\tHG00342\tHG00343\tHG00344\tHG00345\tHG00346\tHG00349\tHG00350\tHG00351\tHG00353\tHG00355\tHG00356\tHG00358\tHG00359\tHG00360\tHG00361\tHG00362\tHG00364\tHG00365\tHG00366\tHG00367\tHG00369\tHG00371\tHG00372\tHG00373\tHG00375\tHG00376\tHG00377\tHG00378\tHG00379\tHG00380\tHG00381\tHG00382\tHG00383\tHG00384\tHG01334\tHG01789\tHG01790\tHG01791\tHG02215\tNA06984\tNA06985\tNA06986\tNA06989\tNA06994\tNA07037\tNA07048\tNA07051\tNA07056\tNA07346\tNA07347\tNA07357\tNA10847\tNA10851\tNA11829\tNA11830\tNA11831\tNA11832\tNA11840\tNA11843\tNA11881\tNA11892\tNA11893\tNA11894\tNA11918\tNA11920\tNA11930\tNA11931\tNA11992\tNA11993\tNA11994\tNA11995\tNA12004\tNA12005\tNA12006\tNA12043\tNA12044\tNA12045\tNA12058\tNA12144\tNA12154\tNA12155\tNA12156\tNA12234\tNA12249\tNA12272\tNA12273\tNA12275\tNA12282\tNA12283\tNA12286\tNA12287\tNA12340\tNA12341\tNA12342\tNA12347\tNA12348\tNA12383\tNA12399\tNA12400\tNA12413\tNA12489\tNA12546\tNA12716\tNA12717\tNA12718\tNA12749\tNA12750\tNA12751\tNA12760\tNA12761\tNA12762\tNA12763\tNA12775\tNA12776\tNA12777\tNA12778\tNA12812\tNA12813\tNA12814\tNA12815\tNA12827\tNA12829\tNA12830\tNA12842\tNA12843\tNA12872\tNA12873\tNA12874\tNA12889\tNA12890\tNA20502\tNA20503\tNA20504\tNA20505\tNA20506\tNA20507\tNA20508\tNA20509\tNA20510\tNA20512\tNA20513\tNA20514\tNA20515\tNA20516\tNA20517\tNA20518\tNA20519\tNA20520\tNA20521\tNA20524\tNA20525\tNA20527\tNA20528\tNA20529\tNA20530\tNA20531\tNA20532\tNA20534\tNA20535\tNA20536\tNA20537\tNA20538\tNA20539\tNA20540\tNA20541\tNA20542\tNA20543\tNA20544\tNA20581\tNA20582\tNA20585\tNA20586\tNA20588\tNA20589\tNA20752\tNA20754\tNA20756\tNA20757\tNA20758\tNA20759\tNA20760\tNA20761\tNA20765\tNA20766\tNA20768\tNA20769\tNA20770\tNA20771\tNA20772\tNA20773\tNA20774\tNA20778\tNA20783\tNA20785\tNA20786\tNA20787\tNA20790\tNA20792\tNA20795\tNA20796\tNA20797\tNA20798\tNA20799\tNA20800\tNA20801\tNA20802\tNA20803\tNA20804\tNA20805\tNA20806\tNA20807\tNA20808\tNA20809\tNA20810\tNA20811\tNA20812\tNA20813\tNA20814\tNA20815\tNA20816\tNA20819\tNA20826\tNA20828" 
#predictionfile_header_eur373="$(head -n 1 ${exprfile_eur} | sed -e 's/,/\t/g')"

predictionfile_header_yri="Gene\tNA18486\tNA18487\tNA18488\tNA18489\tNA18498\tNA18499\tNA18502\tNA18505\tNA18508\tNA18510\tNA18511\tNA18517\tNA18519\tNA18520\tNA18858\tNA18861\tNA18867\tNA18868\tNA18870\tNA18873\tNA18907\tNA18908\tNA18909\tNA18910\tNA18912\tNA18916\tNA18917\tNA18923\tNA18933\tNA18934\tNA19092\tNA19093\tNA19095\tNA19096\tNA19098\tNA19099\tNA19102\tNA19107\tNA19108\tNA19113\tNA19114\tNA19116\tNA19117\tNA19118\tNA19119\tNA19121\tNA19129\tNA19130\tNA19131\tNA19137\tNA19138\tNA19141\tNA19143\tNA19144\tNA19146\tNA19147\tNA19149\tNA19150\tNA19152\tNA19153\tNA19159\tNA19160\tNA19171\tNA19172\tNA19175\tNA19184\tNA19185\tNA19189\tNA19190\tNA19197\tNA19198\tNA19200\tNA19201\tNA19204\tNA19206\tNA19207\tNA19209\tNA19210\tNA19213\tNA19214\tNA19222\tNA19223\tNA19225\tNA19235\tNA19236\tNA19247\tNA19248\tNA19256\tNA19257"
#predictionfile_header_yri="$(head -n 1 ${exprfile_yri89} | sed -e 's/,/\t/g')"

predictionfile_header_eur278="Gene\tHG00096\tHG00097\tHG00099\tHG00100\tHG00101\tHG00102\tHG00103\tHG00104\tHG00105\tHG00106\tHG00108\tHG00109\tHG00110\tHG00111\tHG00112\tHG00114\tHG00115\tHG00116\tHG00117\tHG00118\tHG00119\tHG00120\tHG00121\tHG00122\tHG00123\tHG00124\tHG00125\tHG00126\tHG00127\tHG00128\tHG00129\tHG00130\tHG00131\tHG00132\tHG00133\tHG00134\tHG00135\tHG00136\tHG00137\tHG00138\tHG00139\tHG00141\tHG00142\tHG00143\tHG00145\tHG00146\tHG00148\tHG00149\tHG00150\tHG00151\tHG00152\tHG00154\tHG00155\tHG00156\tHG00157\tHG00158\tHG00159\tHG00160\tHG00231\tHG00232\tHG00233\tHG00234\tHG00235\tHG00236\tHG00238\tHG00239\tHG00240\tHG00242\tHG00243\tHG00244\tHG00245\tHG00246\tHG00247\tHG00249\tHG00250\tHG00251\tHG00252\tHG00253\tHG00255\tHG00256\tHG00257\tHG00258\tHG00259\tHG00260\tHG00261\tHG00262\tHG00263\tHG00264\tHG00265\tHG01334\tHG01789\tHG01790\tHG01791\tHG02215\tNA06984\tNA06985\tNA06986\tNA06989\tNA06994\tNA07037\tNA07048\tNA07051\tNA07056\tNA07346\tNA07347\tNA07357\tNA10847\tNA10851\tNA11829\tNA11830\tNA11831\tNA11832\tNA11840\tNA11843\tNA11881\tNA11892\tNA11893\tNA11894\tNA11918\tNA11920\tNA11930\tNA11931\tNA11992\tNA11993\tNA11994\tNA11995\tNA12004\tNA12005\tNA12006\tNA12043\tNA12044\tNA12045\tNA12058\tNA12144\tNA12154\tNA12155\tNA12156\tNA12234\tNA12249\tNA12272\tNA12273\tNA12275\tNA12282\tNA12283\tNA12286\tNA12287\tNA12340\tNA12341\tNA12342\tNA12347\tNA12348\tNA12383\tNA12399\tNA12400\tNA12413\tNA12489\tNA12546\tNA12716\tNA12717\tNA12718\tNA12749\tNA12750\tNA12751\tNA12760\tNA12761\tNA12762\tNA12763\tNA12775\tNA12776\tNA12777\tNA12778\tNA12812\tNA12813\tNA12814\tNA12815\tNA12827\tNA12829\tNA12830\tNA12842\tNA12843\tNA12872\tNA12873\tNA12874\tNA12889\tNA12890\tNA20502\tNA20503\tNA20504\tNA20505\tNA20506\tNA20507\tNA20508\tNA20509\tNA20510\tNA20512\tNA20513\tNA20514\tNA20515\tNA20516\tNA20517\tNA20518\tNA20519\tNA20520\tNA20521\tNA20524\tNA20525\tNA20527\tNA20528\tNA20529\tNA20530\tNA20531\tNA20532\tNA20534\tNA20535\tNA20536\tNA20537\tNA20538\tNA20539\tNA20540\tNA20541\tNA20542\tNA20543\tNA20544\tNA20581\tNA20582\tNA20585\tNA20586\tNA20588\tNA20589\tNA20752\tNA20754\tNA20756\tNA20757\tNA20758\tNA20759\tNA20760\tNA20761\tNA20765\tNA20766\tNA20768\tNA20769\tNA20770\tNA20771\tNA20772\tNA20773\tNA20774\tNA20778\tNA20783\tNA20785\tNA20786\tNA20787\tNA20790\tNA20792\tNA20795\tNA20796\tNA20797\tNA20798\tNA20799\tNA20800\tNA20801\tNA20802\tNA20803\tNA20804\tNA20805\tNA20806\tNA20807\tNA20808\tNA20809\tNA20810\tNA20811\tNA20812\tNA20813\tNA20814\tNA20815\tNA20816\tNA20819\tNA20826\tNA20828"
#predictionfile_header_eur278="$(head -n 1 ${exprfile_eur278} | sed -e 's/,/\t/g')"

# ==============================================================================================================================
# switches to control which parts of script to run
# set to 0 in order to run
# ==============================================================================================================================
run_compute_weights_eur=1
run_compute_weights_yri=1

run_compute_weights_eur278_toyri=1
run_compute_weights_eur278_tofinn=1
run_compute_weights_eur89_to_yri89=1

run_compute_weights_crosspop=0

# ==============================================================================================================================
# start job scheduling script
# ==============================================================================================================================

# -------------------- #
echo "Start time: $(date)"

# -------------------- #
# make necessary output directories
if [[ ! -d "$logdir" ]]; then mkdir -p $logdir; fi


# -------------------- #
# compute new weights for each gene
# this entails running glmnet on a PLINK RAW genotype dosage file
# and extracting the expression level from the UCSC BED file stored at $exprfile

# how many genes are in this list?
# remember that this file has a header so line count is off by 1 
nGenes=$(wc -l $genelist | cut -f 1 -d " ")
###let "nGenes=nGenes-1"

### uncomment next line for testing and debugging the pipeline
#nGenes=100


# -------------------- #
# first do EUR
if [[ "${run_compute_weights_eur}" -eq "0" ]]; then

    # set variables to population-specific parameters
    # point prediction, weight, lambda, expression, subject files to EUR
    # set alternative population parameters to YRI

    h_rt="23:59:59"
    pop="eur278" 
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

    altpop="yri89"
    altpop_out_lm_file=$out_lm_file_eur2yri
    altpop_out_genelm_file=$out_genelm_file_eur2yri
    altpop_exprfile=$exprfile_yri
    predictionfile_header=$predictionfile_header_eur373
    weightsfile=$weightsfile_eur
    newweightsfile=$newweightsfile_eur
    out_genelm_file=$out_genelm_file_eur
    num_pred_file=$num_pred_file_eur
    nsamples=$nsamples_eur
    out_lm_file=$out_lm_file_eur
    phenofile=$phenofile_eur
    h2file=$h2file_eur
    h2file_null=$h2file_null_eur
    nfolds=$nfolds_eur

    # -------------------- #
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

    # beefy calculations done, don't need as much compute time
    # if possible, then set $h_rt to < 30 min to schedule in QB3 speedy queue 
    h_rt="06:00:00"

    # -------------------- #
    # for each gene, collect weights, lambdas, predictions, etc.
    # merge these results into single file for weights, single file for lambdas, etc.
    qsub -N glmnet.${glmmethod}.collect.weights.${pop} \
         -hold_jid glmnet.${glmmethod}.${pop} \
         -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile,lambdafile=$lambdafile,outdir=$outdir,resultsdir=$resultsdir,resultssubdir=$resultssubdir,tmpdir=$tmpdir,pop=$pop,altpop=$altpop,predictionfile_altpop=$predictionfile_altpop,predictionfile_samepop=$predictionfile_samepop,predictionfile_header=$predictionfile_header,h2file=$h2file,h2file_null=$h2file_null \
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
fi

# -------------------- #
# now do YRI
if [[ "${run_compute_weights_yri}" -eq "0" ]]; then

    # point prediction, weight, lambda, expression, subject files to YRI
    # set alternative population parameters to EUR 
    pop="yri89"
    subjectids_altpop=$subjectids_eur
    h_rt="12:00:00"  # smaller sample size than eur373 --> less compute time necessary
    subjectids=$subjectids_yri
    exprfile=$exprfile_yri
    outdir=$outdir_yri
    resultssubdir=$resultssubdir_yri
    resultsdir=$resultsdir_yri
    subjectids_altpop=$subjectids_eur
    predictionfile=$predictionfile_yri
    lambdafile=$lambdafile_yri
    predictionfile_altpop=$predictionfile_yri2eur
    predictionfile_samepop=$predictionfile_yri2yri

    altpop="eur373"
    altpop_out_lm_file=$out_lm_file_yri2eur
    altpop_out_genelm_file=$out_genelm_file_yri2eur
    altpop_exprfile=$exprfile_eur
    predictionfile_header=$predictionfile_header_yri
    weightsfile=$weightsfile_yri
    newweightsfile=$newweightsfile_yri
    out_genelm_file=$out_genelm_file_yri
    num_pred_file=$num_pred_file_yri
    nsamples=$nsamples_yri
    out_lm_file=$out_lm_file_yri
    phenofile=$phenofile_yri
    h2file=$h2file_yri
    h2file_null=$h2file_null_yri
    nfolds=$nfolds_yri

    # -------------------- #
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
         -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile,lambdafile=$lambdafile,outdir=$outdir,resultsdir=$resultsdir,resultssubdir=$resultssubdir,tmpdir=$tmpdir,pop=$pop,altpop=$altpop,predictionfile_altpop=$predictionfile_altpop,predictionfile_samepop=$predictionfile_samepop,predictionfile_header=$predictionfile_header,h2file=$h2file,h2file_null=$h2file_null \
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
fi


if [[ "${run_compute_weights_eur278_toyri}" -eq "0" ]]; then

    # set variables to population-specific parameters
    # point prediction, weight, lambda, expression, subject files to EUR
    # set alternative population parameters to YRI

    h_rt="23:59:59"
    pop="eur278" 
    subjectids=$subjectids_eur278
    exprfile=$exprfile_eur278
    outdir=$outdir_eur278
    resultssubdir=$resultssubdir_eur278
    resultsdir=$resultsdir_eur278
    predictionfile=$predictionfile_eur278
    predictionfile_header=$predictionfile_header_eur278
    lambdafile=$lambdafile_eur278
    weightsfile=$weightsfile_eur278
    predictionfile_samepop=$predictionfile_eur278toeur278
    nsamples=$nsamples_eur278
    num_pred_file=$num_pred_file_eur278
    out_lm_file=$out_lm_file_eur278
    phenofile=$phenofile_eur278
    newweightsfile=$newweightsfile_eur278
    out_genelm_file=$out_genelm_file_eur278
    h2file=$h2file_eur278
    h2file_null=$h2file_null_eur278
    nfolds=$nfolds_eur278
    
    altpop="yri89"
    predictionfile_altpop=$predictionfile_eur278toyri
    altpop_out_lm_file=$out_lm_file_eur278toyri
    altpop_out_genelm_file=$out_genelm_file_eur278toyri
    subjectids_altpop=$subjectids_yri
    altpop_exprfile=$exprfile_yri

    # -------------------- #
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
         -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile,lambdafile=$lambdafile,outdir=$outdir,resultsdir=$resultsdir,resultssubdir=$resultssubdir,tmpdir=$tmpdir,pop=$pop,altpop=$altpop,predictionfile_altpop=$predictionfile_altpop,predictionfile_samepop=$predictionfile_samepop,predictionfile_header=$predictionfile_header,h2file=$h2file,h2file_null=$h2file_null \
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


fi

if [[ "${run_compute_weights_eur278_tofinn}" -eq "0" ]]; then

    # set variables to population-specific parameters
    # point prediction, weight, lambda, expression, subject files to EUR
    # set alternative population parameters to FINN
    h_rt="23:59:59"
    pop="eur278" 
    subjectids=$subjectids_eur278
    exprfile=$exprfile_eur278
    outdir=$outdir_eur278
    resultssubdir=$resultssubdir_eur278
    resultsdir=$resultsdir_eur278
    predictionfile=$predictionfile_eur278
    predictionfile_header=$predictionfile_header_eur278
    lambdafile=$lambdafile_eur278
    weightsfile=$weightsfile_eur278
    predictionfile_samepop=$predictionfile_eur278toeur278
    nsamples=$nsamples_eur278
    num_pred_file=$num_pred_file_eur278
    out_lm_file=$out_lm_file_eur278
    phenofile=$phenofile_eur278
    newweightsfile=$newweightsfile_eur278
    out_genelm_file=$out_genelm_file_eur278
    h2file=$h2file_eur278
    h2file_null=$h2file_null_eur278
    nfolds=$nfolds_eur278

    
    altpop="fin95"
    altpop_exprfile=$exprfile_fin
    subjectids_altpop=$subjectids_fin
    predictionfile_altpop=$predictionfile_eur278tofin
    altpop_out_lm_file=$out_lm_file_eur278tofin
    altpop_out_genelm_file=$out_genelm_file_eur278tofin

    # -------------------- #
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
         -v weightsfile=$weightsfile,glmnetdir=$glmnetdir,glmmethod=$glmmethod,logdir=$logdir,predictionfile=$predictionfile,lambdafile=$lambdafile,outdir=$outdir,resultsdir=$resultsdir,resultssubdir=$resultssubdir,tmpdir=$tmpdir,pop=$pop,altpop=$altpop,predictionfile_altpop=$predictionfile_altpop,predictionfile_samepop=$predictionfile_samepop,predictionfile_header=$predictionfile_header,h2file=$h2file,h2file_null=$h2file_null \
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


fi



# this scenario performs subsampling of EUR373, pulling 89 people per subsample
if [[ "${run_compute_weights_eur89_to_yri89}" -eq "0" ]]; then


    # ensure that eur89_dir exists
    mkdir -p $eur89_dir

    # perform subsampling once, prior to constructing predicion models
    # arg order:
    # 1) number of resamplings to perform
    # 2) file path to EUR373 expression data
    # 3) directory path to EUR89 output directory
    #
    # output, for each subsample $i:
    # 1) a TXT file of RNA data, ngenes x 89, in ${eur89_dir}/geuvadis.eur89_${i}.RPKM.invnorm.txt
    # 2) a PLINK PHENO file of RNA data, 89 x nGenes, in ${eur89_dir}/geuvadis.eur89_${i}.RPKM.invnorm.pheno
    # 3) a PLINK FID/IID file of subject IDs, 89 x 2, in ${eur89_dir}/geuvadis.eur89_${i}.subjectids.txt
    $Rscript $R_subsample_eur $nresample $exprfile_eur $eur89_dir 

    # loop over all subsamples
    # unlike previous setup, will need to manually write the variables particular to each eur89
    for i in $(seq 1 $nresample); do

        # $popi is merely shorthand for current subsample
        popi="${pop}_${i}"

        # all of these filepaths depend on $i through $popi
        subjectids_eur89="${eur89_dir}/geuvadis.${popi}.sampleids.txt"
        exprfile_eur89="${eur89_dir}/geuvadis.${popi}.RPKM.invnorm.txt"
        outdir_eur89="/scratch/kkeys/${glmmethod}/genes/eur89/${popi}"
        resultsdir_eur89="${glmnetdir}/${glmmethod}/eur89/${popi}"
        resultssubdir_eur89="${resultsdir_eur89}/results"
        predictionfile_eur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_predictions.txt"
        predictionfile_header_eur89="$(head -n 1 ${exprfile_eur89} | sed -e 's/,/\t/g')"
        lambdafile_eur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_lambdas.txt"
        weightsfile_eur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_weights.txt"
        predictionfile_eur89toeur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_predictinto_${popi}.txt"
        num_pred_file_eur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_numpred.txt"
        out_lm_file_eur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_lm_predvmeas_results.txt"
        phenofile_eur89="${eur89_dir}/geuvadis.${popi}.RPKM.invnorm.pheno"
        newweightsfile_eur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_weights_noNA.txt"
        out_genelm_file_eur89="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_genelm_predvmeas_results.txt"
        h2file_eur89="${resultsdir_eur89}/geuvadis_h2_${popi}.txt"
        h2file_null_eur89="${resultsdir_eur89}/geuvadis_h2_null_${popi}.txt"

        predictionfile_eur89toyri="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_predictinto_yri89.txt"
        out_lm_file_eur89toyri="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_predictinto_yri89_lm_predvmeas_results.txt"
        out_genelm_file_eur89toyri="${resultsdir_eur89}/geuvadis_${glmmethod}_${popi}_predictinto_yri89_lm_predvmeas_results.txt"

        # ensure that all output directories exist
        mkdir -p $outdir_eur89
        mkdir -p $resultssubdir_eur89

        # this section of code seems redundant, but really it's just lazy
        # easier to include this than to manually change all variables in QSUB calls below
        subjectids=$subjectids_eur89
        exprfile=$exprfile_eur89
        outdir=$outdir_eur89
        resultssubdir=$resultssubdir_eur89
        resultsdir=$resultsdir_eur89
        predictionfile=$predictionfile_eur89
        predictionfile_header=$predictionfile_header_eur89
        lambdafile=$lambdafile_eur89
        weightsfile=$weightsfile_eur89
        predictionfile_samepop=$predictionfile_eur89toeur89
        nsamples=$nsamples_eur89
        num_pred_file=$num_pred_file_eur89
        out_lm_file=$out_lm_file_eur89
        phenofile=$phenofile_eur89
        newweightsfile=$newweightsfile_eur89
        out_genelm_file=$out_genelm_file_eur89
        h2file=$h2file_eur89
        h2file_null=$h2file_null_eur89
        nfolds=$nfolds_eur89

        # altpop is always YRI 
        altpop="yri89"
        predictionfile_altpop=$predictionfile_eur89toyri
        altpop_out_lm_file=$out_lm_file_eur89toyri
        altpop_out_genelm_file=$out_genelm_file_eur89toyri
        subjectids_altpop=$subjectids_yri
        altpop_exprfile=$exprfile_yri

        # -------------------- #
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

    done

fi

# ==============================================================================================================================
# cross-population training
# ==============================================================================================================================
#
# in this section we will train from CEU, GBR, TSI, FIN, YRI into all other populations
# for example, we train in CEU and predict into other four populations
#
# key is that we downsample to smallest population. sample sizes:
# -- YRI = 89
# -- CEU = 92
# -- FIN = 95
# -- GBR = 96
# -- TSI = 93
#
# thus, we must downsample to 89 samples in each training population 

if [[ "${run_compute_weights_crosspop}" -eq "0" ]]; then

    # arrays of pops + their respective sample sizes
    pops=("ceu" "tsi" "gbr" "fin" "yri")
    popsizes=("92" "93" "96" "95" "89")

    # set runtime to relatively high value
    h_rt="23:59:59"
    #h_rt="00:29:59"

    # loop over pops
    for i in $(seq 0 4); do 

        # set variables for current pop 
        # train in "pop", test in "notpop" 
        pop="${pops[$i]}"  ## note absence of sample size; will add subsample size 89 to this later 
        notpop="not${pop}"
        popsize="${popsizes[$i]}"
        altpop="${notpop}"
        subjectids="${geuvadisdir}/geuvadis.${pop}${popsize}.sampleids.txt"

        # first subsample the population in question
        $Rscript $R_subsample_pop $pop $exprfile_eur $exprfile_yri $subjectids $crosspop_dir

        # rename population with correct subsampled population size
        pop="${pop}89"

        outdir="/scratch/kkeys/${glmmethod}/genes/${pop}"
        resultsdir="${resultsdir_crosspop}/${pop}"
        resultssubdir="${resultsdir}/results"

        mkdir -p $outdir
        mkdir -p $resultsdir
        mkdir -p $resultssubdir

        subjectids="${crosspop_dir}/geuvadis.${pop}.sampleids.txt"
        exprfile="${crosspop_dir}/geuvadis.${pop}.RPKM.invnorm.txt"
        predictionfile="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictions.txt"
        lambdafile="${resultsdir}/geuvadis_${glmmethod}_${pop}_lambdas.txt"
        weightsfile="${resultsdir}/geuvadis_${glmmethod}_${pop}_weights.txt"
        predictionfile_samepop="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${pop}.txt"
        num_pred_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_numpred.txt"
        nsamples="89"
        out_lm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_lm_predvmeas_results.txt"
        newweightsfile="${resultsdir}/geuvadis_${glmmethod}_${pop}_weights_noNA.txt"
        out_genelm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_genelm_predvmeas_results.txt"
        h2file="${resultsdir}/geuvadis_h2_${pop}.txt"
        h2file_null="${resultsdir}/geuvadis_h2_null_${pop}.txt"
        nfolds=$nfolds_yri
        phenofile="${crosspop_dir}/geuvadis.${pop}.RPKM.invnorm.pheno"

        altpop_exprfile="${crosspop_dir}/geuvadis.${notpop}.RPKM.invnorm.txt"
        subjectids_altpop="${crosspop_dir}/geuvadis.${notpop}.sampleids.txt"
        predictionfile_altpop="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${notpop}.txt"
        altpop_out_lm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${notpop}_lm_predvmeas_results.txt"
        altpop_out_genelm_file="${resultsdir}/geuvadis_${glmmethod}_${pop}_predictinto_${notpop}_genelm_predvmeas_results.txt"

### TODO: Ensure that $predictionfile_header causes no problems!!!!
        predictionfile_header="$(head -n 1 ${exprfile} | sed -e 's/,/\t/g')"


        source geuvadis_qsub_jobs.sh
    done
fi

# end script
echo "End time: $(date)"
