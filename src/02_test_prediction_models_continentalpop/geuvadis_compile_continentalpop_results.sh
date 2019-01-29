#!/usr/bin/env bash

Rscript=${Rscript}

R_compile_results_continentalpop="${R_compile_results_continentalpop}"

eur373_to_eur373_path=""
eur373_to_afr_path=""
eur278_to_eur278_path=""
eur278_to_fin_path=""
eur278_to_afr_path=""
afr_to_eur373_path=""
afr_to_eur278_path=""
afr_to_fin_path=""
afr_to_afr_path=""
output_results=""
output_r2=""

$RSCRIPT $R_compile_results \
    --EUR373-to-EUR373 ${eur373_to_eur373_path} \
    --EUR373-to-AFR ${eur373_to_afr_path} \
    --EUR278-to-EUR278 ${eur278_to_eur278_path} \
    --EUR278-to-FIN ${eur278_to_fin_path} \
    --EUR278-to-AFR ${eur278_to_afr_path} \
    --AFR-to-EUR373 ${afr_to_eur373_path} \
    --AFR-to-EUR278 ${afr_to_eur278_path} \
    --AFR-to-FIN ${afr_to_fin_path} \
    --AFR-to-AFR ${afr_to_afr_path} \
    --output-results ${output_results} \
    --output-r2 ${output_r2}
