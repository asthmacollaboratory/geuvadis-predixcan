#!/usr/bin/env bash
# ==========================================================================================
# coded by Kevin L. Keys (2018)
#
# This script parses GEUVADIS expression and genotype data. It produces the following files:
#
#   (1) headered and labeled matrices of expression levels, plus their transposes,
#       and lists of sample IDs, for each of EUR and YRI populations
#   (2) BGZIP'd and tabix-indexed VCFs of the gEUVADIS genotype VCFs
#   (3) a filtered, postprocessed trio of PLINK BED/BIM/FAM that covers all chromosomes 
#
# ==========================================================================================

# ==========================================================================================
# BASH script settings
# ==========================================================================================
set -e  ## script will exit on error
set -u  ## script will exit if it sees an uninitialized variable

# ==========================================================================================
# directories 
# ==========================================================================================
thisdir="$(dirname $(readlink -f $0))"

analysisdir="${thisdir}/analysis"
datadir="${analysisdir}/data"
geuvadis_genodir="${datadir}/genotypes"
geuvadis_rnaseqdir="${datadir}/rnaseq"
vcfdir="${geuvadis_genodir}/vcf"
plinkdir="${geuvadis_genodir}/plink"


# ==========================================================================================
# binaries
# ==========================================================================================
# these are paths to static executables
# note that these are GUESSED with "whereis". this is a *fragile* way to set these variables!
# "whereis" can return an empty result. the guess also uses the first result, which may not be desired. 
# a better solution, if it is available, is to point these variables directly to desired binaries
R=$(whereis R | awk '{print $2}')
RSCRIPT=$(whereis Rscript | awk '{print $2}')
PLINK=$(whereis plink | awk '{print $2}')  ## <-- note: by default, doesn't use PLINK2
PYTHON=$(whereis python3 | awk '{print $2}')
GCTA=$(whereis gcta64 | awk '{print $2}')
VCFFILTER=$(whereis vcffilter| awk '{print $2}')
BCFTOOLS=$(whereis bcftools| awk '{print $2}')

# ==========================================================================================
# external scripts 
# ==========================================================================================
R_parse_geu_expr="${thisdir}/parse_geuvadis_expressions.R"
R_remap_snpids="${thisdir}/convert_geuvadis_snpids_to_rsids.R" 


# ==========================================================================================
# file paths 
# ==========================================================================================
# static file paths 
vcfdatapfx="GEUVADIS"
vcfdatasfx="PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes"


# script file paths
all_rnaseq_data_file="${geuvadis_rnaseqdir}/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt"
sample_ids_file="${geuvadis_genodir}/sample_ids/gevuadis.465.sample.ids.txt"
ryan_file="${geuvadis_rnaseqdir}/expressions.genes.yri.matrix.eqtl.txt"
eur_out_file="${geuvadis_rnaseqdir}/geuvadis.eur373.RPKM.invnorm.txt"
yri_out_file="${geuvadis_rnaseqdir}/geuvadis.yri89.RPKM.invnorm.txt"
teur_out_file="${geuvadis_rnaseqdir}/geuvadis.eur373.RPKM.invnom.transposed.txt"
teur_out_file="${geuvadis_rnaseqdir}/geuvadis.eur373.RPKM.invnom.transposed.txt"
eur_ids_file="${geuvadis_rnaseqdir}/geuvadis.eur373.ids.txt"
yri_ids_file="${geuvadis_rnaseqdir}/geuvadis.yri373.ids.txt"
mergelist="${geuvadis_genodir}/plink/geuvadis_mergelist_for_plink.txt"
id_map_file="${geuvadis_genodir}/snp_ids/Phase1.Geuvadis_dbSnp137_idconvert.txt"


# ==========================================================================================
# script variables, URLs 
# ==========================================================================================
nthreads=12     # this applies for each chromosome
mind="0.05"     # threshold for inclusion of missingness per sample
geno="0.05"     # threshold for inclusion of missingness per genotype
maf="0.01"      # threshold for minimum permissible minor allele frequency
hwe="0.00001"   # purge markers whose p-value for deviation from Hardy-Weinberg equilibrium falls below this number
rsq="0.3"       # threshold for minimum permissible imputation R-squared 

# ==========================================================================================
# script switches 
# ==========================================================================================
# toggle various parts of the pipeline here 
# set appropriate variable to 0 to run
# any other value turns that analysis off
run_parse_expr_data=1
run_parse_genotypes=0

# ==========================================================================================
# executable code 
# ==========================================================================================

# first make directories, if they do not already exist
mkdir -p ${analysisdir}
mkdir -p ${datadir}
mkdir -p ${geuvadis_genodir}
mkdir -p ${geuvadis_rnaseqdir}
mkdir -p ${vcfdir}
mkdir -p ${plinkdir}

### --- download genotype data --- ###
for chr in $(seq 1 22); do
    current_vcf_url="https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GEUVADIS.chr${chr}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
    wget ${current_vcf_url} -P ${vcfdir} 
done

### --- download expression data --- ###
rnaseq_url="https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
wget ${rnaseq_url} -P ${rnaseqdir}

### --- parse expression data --- ### 
if [[ "${run_parse_expr_data}" -eq "0" ]]; then

    # everything is in one script
    $RSCRIPT $R_parse_geu_expr ${all_rnaseq_data_file} ${sample_ids_file} ${ryan_file} ${eur_out_file} ${yri_out_file} ${teur_out_file} ${tyri_out_file} ${eur_ids_file} ${yri_ids_file}

fi


### --- parse genotype data --- ### 
if [[ "${run_parse_genotypes}" -eq "0" ]]; then


    # ensure that the VCFs are all bgzip'd and tabix'd
    for file in $(ls ${genovadis_genodir}/vcf/${vcfdatapfx}.*.vcf.gz); do
        mv -f $file backup.vcf;
        $BCFTOOLS view --output-type z --threads $nthreads -o $file backup.vcf;
        ~/bin/htsfile $file;
        $BCFTOOLS index --tbi $file --threads $nthreads;
    done

    # set up a mergelist for PLINK
    rm -f $mergelist
    touch $mergelist

    # loop over each chromosome
    for chr in $(seq 1 22); do
        
        # file names that depend on chr
        vcfin="${vcfdir}/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.vcf.gz"
        vcfout="${vcfdir}/${vcfdatapfx}.chr${chr}.${vcfdatasfx}_rsq_${rsq}.vcf.gz"
        plinkout="${plinkdir}/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.rsq_${rsq}_maf_${maf}_hwe_${hwe}_geno_${geno}"

        # save current PLINK file to merge
        echo $plinkout >> $mergelist

        # now filter by imputation R2
        $BCFTOOLS filter --include "INFO/RSQ > ${rsq} & INFO/VT='SNP'" --output-type z --output ${vcfout} --threads ${nthreads} ${vcfin}

        # make sure to reindex with tabix protocol
        $BCFTOOLS index  --threads $nthreads --tbi ${vcfout}
    
        # convert VCF to PLINK format
        # apply filters along the way
        $PLINK \
            --vcf $vcfout \
            --keep-allele-order \
            --vcf-idspace-to ":" \
            --biallelic-only strict list \
            --double-id \
            --hwe $hwe midp \
            --maf $maf \
            --geno $geno \
            --make-bed \
            --out $plinkout \
            --threads $nthreads

        # backup the BIM file
        cp ${plinkout}.bim ${plinkout}.bim.backup

        # replace gEUVADIS SNP IDs with rsIDs 
        # this will read from the *backup* BIM file and clobber the *original* BIM 
        $RSCRIPT $R_remap_snpids ${plinkout}.bim.backup ${id_map_file} ${plinkout}.bim

    done

    # new PLINK output prefix 
    plinkout="${plinkdir}/${vcfdatapfx}.ALLCHR.${vcfdatasfx}.rsq_0.8_maf_${maf}_hwe_${hwe}_geno_${geno}"

    # now merge all BEDs into a single PLINK set
    $PLINK \
        --merge-list $mergelist \
        --make-bed \
        --out $plinkout \
        --threads $nthreads
fi

exit 0

### old code
#        # replace gEUVADIS SNP IDs with rsIDs 
#        # this will read from the *backup* BIM file and clobber the *original* BIM 
#        while read -r -a line; do
#
#            # parse columns of BIM
#            chr="${line[0]}"
#            id="${line[1]}"
#            gendist="${line[2]}"
#            pos="${line[3]}"
#            A1="${line[4]}"
#            A2="${line[5]}"
#            
#            # find rsID in the gEUVADIS SNP ID conversion file
#            rsid=$(fgrep "${id}" ../snp_ids/Phase1.Geuvadis_dbSnp137_idconvert.txt | awk '{print $1}')
#
#            # write 
#            echo -e "${chr}\t${rsid}\t${gendist}\t${pos}\t${A1}\t${A2}"
#        done < ${plinkout}.bim.backup > ${plinkout}.bim 

#
#
#cat << done-catting > sub_parse_geuvadis_data_chr${chr}.sh
#        
#        # file names that depend on chr
#        vcfin="${geuvadis_genodir}/vcf/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.vcf.gz"
#        vcfout="${geuvadis_genodir}/vcf/${vcfdatapfx}.chr${chr}.${vcfdatasfx}_rsq_${rsq}.vcf.gz"
#        plinkout="${geuvadis_genodir}/plink/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.rsq_0.8_maf_${maf}_hwe_${hwe}_geno_${geno}"
#
#        # save current PLINK file to merge
#        echo $plinkout >> $mergelist
#
#        # now filter by imputation R2
#        # make sure to reindex with tabix protocol
#        bcftools filter --include "INFO/RSQ > ${rsq} & INFO/VT='SNP'" --output-type z --output ${vcfout} --threads ${nthreads} ${vcfin}
#        bcftools index  --threads $nthreads --tbi ${vcfout}
#    
#        # convert VCF to PLINK format
#        # apply filters along the way
#        $PLINK \\
#            --vcf $vcfout \\
#            --keep-allele-order \\
#            --vcf-idspace-to ":" \\
#            --biallelic-only strict list \\
#            --double-id \\
#            --hwe $hwe midp \\
#            --maf $maf \\
#            --geno $geno \\
#            --make-bed \\
#            --out $plinkout \\
#            --threads $nthreads
#
#        # backup the BIM file
#        cp ${plinkout}.bim ${plinkout}.bim.backup
#
#        
#done-catting
