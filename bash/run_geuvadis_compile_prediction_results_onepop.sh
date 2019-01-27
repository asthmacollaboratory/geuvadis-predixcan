#!/usr/bin/env bash

crosspop_dir="$HOME/gala_sage/rnaseq/geuvadis/crosspop"
prediction_dir="$HOME/gala_sage/rnaseq/glmnet/elasticnet/crosspop"
outdir="$HOME/tmp"
#outdir=${crosspop_dir}
output_pfx="${outdir}/geuvadis.crosspop.results"


pops=("ceu" "tsi" "gbr" "yri" "fin")
altpops=("notceu" "nottsi" "notgbr" "notyri" "notfin")

popkey="$HOME/gala_sage/rnaseq/geuvadis/crosspop/geuvadis.crosspop.sample.pop.key.txt"

RSCRIPT="$HOME/bin/Rscript"
#R_compile_results="$HOME/gala_sage/code/geuvadis_compile_prediction_results_onepop.R"
R_compile_results="/netapp/home/kkeys/Git/geuvadis-predixcan/R/geuvadis_compile_prediction_results_onepop.R"

#for i in $(seq 0 4); do
for i in $(seq 4 4); do

    pop="${pops[$i]}"
    altpop="${altpops[$i]}"

    rna_pop="${crosspop_dir}/geuvadis.${pop}89.RPKM.invnorm.txt"
    rna_altpop="${crosspop_dir}/geuvadis.${altpop}.RPKM.invnorm.txt"
    prediction_popdir="${prediction_dir}/${pop}89"
    predictions_pop="${prediction_popdir}/geuvadis_elasticnet_${pop}89_predictions.txt"
    predictions_altpop="${prediction_popdir}/geuvadis_elasticnet_${pop}89_predictinto_${altpop}.txt"

    nohup $RSCRIPT $R_compile_results \
        --pop ${pop} \
        --altpop ${altpop} \
        --rna-pop ${rna_pop} \
        --rna-altpop ${rna_altpop} \
        --predictions-pop ${predictions_pop} \
        --predictions-altpop ${predictions_altpop} \
        --sample-pop-key ${popkey} \
        --output-prefix ${output_pfx} \
        > ${outdir}/nohup.geuvadis.results.${pop}.out \
        2> ${outdir}/nohup.geuvadis.results.${pop}.err &
done
