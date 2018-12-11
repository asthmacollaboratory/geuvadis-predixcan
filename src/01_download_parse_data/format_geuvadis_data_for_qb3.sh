#!/usr/bin/env bash

# fix a HOME path
# this enables non-Kevins to run this script from their personal $HOME
MYHOME="/media/BurchardRaid01/LabShare/Home/kkeys"

# static executable paths
R="${MYHOME}/bin/Rmkl"
Rscript="${MYHOME}/bin/Rscriptmkl"
PLINK="${MYHOME}/bin/plink"
PYTHON="/usr/local/bin/python3"
GCTA="${MYHOME}/bin/gcta64"
VCFFILTER="${MYHOME}/bin/vcffilter"

# static file paths 
vcfdir="${MYHOME}/gala_sage/rnaseq/geuvadis/genotypes"
vcfdatapfx="GEUVADIS"
vcfdatasfx="PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes"
#lab2nwd="/media/BurchardRaid01/LabShare/Home/cmak/wgs/joint/FINAL.2015-10-30/data/lab2nwd.txt"
#phenofile="${MYHOME}/gala_sage/phenotypes/gala_sage_phenotypes/sage2_clean2015_03_10_de_ident.csv"

# script directories
rnaseqdir="${MYHOME}/gala_sage/rnaseq/geuvadis/rnaseq"
datadir="${rnaseqdir}/data"
#codedir="${rnaseqdir}/code"
codedir="${MYHOME}/gala_sage/rnaseq/code"
resultsdir="${rnaseqdir}/results"
genodir="${MYHOME}/gala_sage/genotypes/mergedLAT-LATP/SAGE_merge"
fastqtldir="${MYHOME}/Git/fastqtl"
galasagedir="${MYHOME}/gala_sage"
rnaseqdir="${galasagedir}/rnaseq"
geuvadisdir="${rnaseqdir}/geuvadis"
geuvadis_genodir="${geuvadisdir}/genotypes"
geuvadis_rnaseqdir="${geuvadisdir}/rnaseq"

# script file paths
all_rnaseq_data_file="${geuvadis_rnaseqdir}/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt"
sample_ids_file="${geuvadis_genodir}/sample_ids/gevuadis.465.sample.ids.txt"
snp_ids_file="${geuvadis_genodir}/snp_ids/Phase1.Geuvadis_dbSnp137_idconvert.txt"
ryan_file="${geuvadis_rnaseqdir}/expressions.genes.yri.matrix.eqtl.txt"
eur_out_file="${geuvadis_rnaseqdir}/geuvadis.eur373.RPKM.invnorm.txt"
yri_out_file="${geuvadis_rnaseqdir}/geuvadis.yri89.RPKM.invnorm.txt"
teur_out_file="${geuvadis_rnaseqdir}/geuvadis.eur373.RPKM.invnom.transposed.txt"
teur_out_file="${geuvadis_rnaseqdir}/geuvadis.eur373.RPKM.invnom.transposed.txt"
eur_ids_file="${geuvadis_genodir}/sample_ids/geuvadis.eur373.sampleids.txt"
yri_ids_file="${geuvadis_genodir}/sample_ids/geuvadis.yri89.sampleids.txt"
mergelist="${geuvadis_genodir}/plink/geuvadis_mergelist_for_plink.txt"
allchrfile="${geuvadis_genodir}/plink/${vcfdatapfx}.ALLCHR.${vcfdatasfx}.rsq_0.8_maf_${maf}_hwe_${hwe}_geno_${geno}"
qb3_geuvadis_genodir="/netapp/home/kkeys/gala_sage/genotypes/gEUVADIS"
qb3_geuvadis_rnaseqdir="/netapp/home/kkeys/gala_sage/rnaseq/geuvadis"
genechrpos_file="${geuvadis_rnaseqdir}/human_ens_GRCh37_genechrpos.txt"
master_genechrpos_file="${MYHOME}/gala_sage/rnaseq/data/human_ens_GRCh37_annot.extended.txt"
genechrpos_500kbfile="${geuvadis_rnaseqdir}/human_ens_GRCh37_genechrpos_plusmin500kb.txt"

# script names
R_parse_geu_expr="${codedir}/parse_geuvadis_expressions.R"
R_remap_snpids="${codedir}/convert_geuvadis_snpids_to_rsids.R"

# other script variables
nthreads=12     # this applies for each chromosome
mind="0.05"     # threshold for inclusion of missingness per sample
geno="0.05"     # threshold for inclusion of missingness per genotype
maf="0.01"      # threshold for minimum permissible minor allele frequency
hwe="0.00001"   # purge markers whose p-value for deviation from Hardy-Weinberg equilibrium falls below this number
rsq="0.3"       # threshold for minimum permissible imputation R-squared 

if [[ 1 -lt 0 ]]; then

for chr in $(seq 1 22);
do

    # file names that depend on chr
    vcfin="${geuvadis_genodir}/vcf/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.vcf.gz"
    vcfout="${geuvadis_genodir}/vcf/${vcfdatapfx}.chr${chr}.${vcfdatasfx}_rsq_${rsq}.vcf.gz"
    plinkout="${geuvadis_genodir}/plink/${vcfdatapfx}.chr${chr}.${vcfdatasfx}.rsq_0.8_maf_${maf}_hwe_${hwe}_geno_${geno}"

    # save current PLINK file to merge
    #echo $plinkout >> $mergelist

    # now filter by imputation R2
    # make sure to reindex with tabix protocol
    #bcftools filter --include "INFO/RSQ > ${rsq} & INFO/VT='SNP'" --output-type z --output ${vcfout} --threads ${nthreads} ${vcfin}
    #bcftools index  --threads $nthreads --tbi ${vcfout}

    # convert VCF to PLINK format
    # apply filters along the way
    $PLINK \
    --vcf $vcfout \
    --keep-allele-order \
    --vcf-idspace-to ":" \
    --biallelic-only strict list\
    --double-id \
    --make-bed \
    --out $plinkout \
    --threads $nthreads
    ### 10 April 2018: this VCF contains both EUR and YRI, so we should not apply MAF/HWE/missingness filters yet 
    ### since this includes 1KG samples, then we may not need to QC the genotypes at all
    #--hwe $hwe midp \
    #--maf $maf \
    #--geno $geno \

    # backup the BIM file
    rm -f ${plinkout}.bim.backup
    cp -f ${plinkout}.bim ${plinkout}.bim.backup

    # replace gEUVADIS SNP IDs with rsIDs 
    # this will read from the *backup* BIM file and clobber the *original* BIM 
    $Rscript $R_remap_snpids ${plinkout}.bim.backup $snp_ids_file ${plinkout}.bim
done

# merge PLINK files into one big genomic file
$PLINK \
    --merge-list $mergelist \
    --make-bed \
    --out $allchrfile 

scp -i $qb3key ${allchrfile}.* $qb3:$qb3_geuvadis_genodir

fi

### make files for sample IDs
IFS=',' read -r -a eurnames <<< $(head -n 1 $eur_out_file | sed -e "s/Gene,//")
for i in "${eurnames[@]}"; do
    echo "$i $i"
done > $eur_ids_file

IFS=',' read -r -a yrinames <<< $(head -n 1 $yri_out_file | sed -e "s/Gene,//")
for i in "${yrinames[@]}"; do
    echo "$i $i"
done > $yri_ids_file

scp -i $qb3key $eur_ids_file $qb3:$qb3_geuvadis_rnaseqdir
scp -i $qb3key $yri_ids_file $qb3:$qb3_geuvadis_rnaseqdir



### parse gene-chr-start-end info
### then add 500kb up/downstream to gene 
cut -f 1,5,6,7 $master_genechrpos_file > $genechrpos_file 
Rscript \
    -e "library(data.table)" \
    -e "x = fread(\"$genechrpos_file\")" \
    -e "names(x)[1:2] = c(\"Gene\", \"chr\")" \
    -e "chr1to22 = x\$chr %in% 1:22" \
    -e "x = x[chr1to22,]" \
    -e "x\$start_pos = x\$Start_Position - 500000" \
    -e "x\$start_pos[x\$start_pos < 0] = 1" \
    -e "x\$end_pos = x\$End_Position + 500000" \
    -e "y = subset(x, select = c(\"Gene\", \"chr\", \"start_pos\", \"end_pos\"))" \
    -e "fwrite(y, file = \"$genechrpos_500kbfile\", sep = \" \", row.names = F, col.names = F, quote = F)"

scp -i $qb3key $genechrpos_500kbfile $qb3:$qb3_geuvadis_rnaseqdir

### send normalized RNA measures
### we will keep EUR and YRI separate
scp -i $qb3key $eur_out_file $qb3:$qb3_geuvadis_rnaseqdir
scp -i $qb3key $yri_out_file $qb3:$qb3_geuvadis_rnaseqdir
