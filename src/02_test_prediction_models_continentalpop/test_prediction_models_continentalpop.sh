#!/usr/bin/env bash
# ==========================================================================================
# copyright Asthma Collboratory (2018)
# coded by Kevin L. Keys
#
# This script computes prediction weights from GEUVADIS transcriptome data.
# The prediction models are estimated in accordance with PrediXcan protocols;
# for reference, see:
# -- Gamazon et al. (2015) Nature Genetics
# -- https://github.com/hakyimlab/PrediXcan/tree/master/Paper-Scripts
# 
# The GEUVADIS populations used here are:
# 1. EUR373, which includes all 373 unrelated Europeans from GEUVADIS
# 2. EUR278, which includes all 278 non-Finnish Europeans
# 3. FIN, which includes all 95 Finns from GEUVADIS
# 4. YRI, which includes all 89 Yoruba from GEUVADIS
#
# This script schedules four train-test scenarios:
# -- EUR373 to YRI
# -- YRI to EUR373
# -- EUR278 to YRI
# -- EUR278 to FIN
#
# Call:
#
# ./test_prediction_models_continentalpop.sh $ALPHA
#
# where
# -- $ALPHA = "0.0", "0.5", or "1.0", used by glmnet to determine the training algorithm
# ==========================================================================================


# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable


# ==========================================================================================
# parse command line arguments
# ==========================================================================================

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


# ==========================================================================================
# script variables
# ==========================================================================================

thisdir="$(dirname $(readlink -f $0))"

BASH_schedule_eur373_to_yri="${thisdir}/schedule_eur373_to_yri.sh"
BASH_schedule_yri_to_eur373="${thisdir}/schedule_yri_to_eur373.sh"
BASH_schedule_eur278_to_fin="${thisdir}/schedule_eur278_to_fin.sh"
BASH_schedule_eur278_to_yri="${thisdir}/schedule_eur278_to_yri.sh"
BASH_define_variables="${thisdir}/../common/geuvadis_variables.sh"


# ==============================================================================================================================
# source variables 
# ==============================================================================================================================

# must first grab all of the variables for GEUVADIS analysis 
source ${BASH_define_variables}


# ==========================================================================================
# start job scheduling script
# ==========================================================================================

# -------------------- #
# EUR373 to YRI 

source ${BASH_schedule_eur373_to_yri}

# -------------------- #
# YRI to EUR373 

source ${BASH_schedule_yri_to_eur373}

# -------------------- #
# EUR278 to YRI

source ${BASH_schedule_eur278_to_yri}

# -------------------- #
# EUR278 to FIN

source ${BASH_schedule_eur278_to_fin}
