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
# ./test_prediction_models_continentalpop.sh
# ==========================================================================================


# ==========================================================================================
# BASH script settings
# ==========================================================================================

set -o errexit # set -e, script will exit on error
set -o nounset # set -u, script will exit if it sees an uninitialized variable
#set -o xtrace  # set -x, script will track which command is currently running

# ==========================================================================================
# script variables
# ==========================================================================================

thisdir="$(dirname $(readlink -f $0))"

BASH_schedule_eur373_to_yri="${thisdir}/schedule_eur373_to_yri.sh"
BASH_schedule_yri_to_eur373="${thisdir}/schedule_yri_to_eur373.sh"
BASH_schedule_eur278_to_fin="${thisdir}/schedule_eur278_to_fin.sh"
BASH_schedule_eur278_to_yri="${thisdir}/schedule_eur278_to_yri.sh"
BASH_define_variables="${thisdir}/../common/geuvadis_variables.sh"


# ==========================================================================================
# source variables
# ==========================================================================================

# must first grab all of the variables for GEUVADIS analysis
source ${BASH_define_variables}


# ==========================================================================================
# start job scheduling script
# ==========================================================================================

# -------------------- #
# EUR373 to YRI

source ${BASH_schedule_eur373_to_yri}
#
## -------------------- #
## YRI to EUR373
#
#source ${BASH_schedule_yri_to_eur373}
#
## -------------------- #
## EUR278 to YRI
#
#source ${BASH_schedule_eur278_to_yri}
#
## -------------------- #
## EUR278 to FIN
#
#source ${BASH_schedule_eur278_to_fin}
