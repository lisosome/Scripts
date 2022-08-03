#!/bin/env bash

# script to perform sleep/lipids GWAS
# will run {M1,M2} x {MALE,FEMALE,COMBINED} x {LTST,STST} in single run
# Assuming run by chromosome so

# set to chromosome. Used in output file name to prevent output being clobbered
chr=$1
trt=$2
out_dir=$3
# path to MMAP
mmap=/home/nardone/software/mmap/mmap_gli_training/programs/mmap.2022_02_08_14_53.intel

# path to MMAP binary genotype file
geno=/netapp05/analisi_nardone/CHARGE_DeprBP/binary_files/genotype_binary/INGI-FVG.vcf_converted.1.r2.0.4.mac.5.0.bit.12.probt.0.99.b.1001.chr${chr}.bin

# path to MMAP phenotype file
pheno=/netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/2.INGI-FVG.pheno.1927.sex_combined.DBP.csv

# set trait name
trait=${trt}

#set covariates
# !! DO NOT ADD sex
covariates=" PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE2 "

# set the strata variable sex
# !! sex MUST BE 0/1
#strata_str="  "
strata_str=" --strata SEX "

# set ExC terms
# these may be added to the phenotype file directly but
# MMAP will not use them
# !! DO NOT ADD sex
ExC_str=" "
ExC_str=" --ExC_covariates AGE AGE2  "

# can drop ExC covariates in general but not for GLI
exstrata_str=" "
exstrata_str=" --include_ExStrata "

exposures=" dDEPR qDEPR "
#
exposure_str=" "
exposure_str=" --exposures ${exposures} "

# if set will exit before GWAS to be able to examine model fit
#null_model_str=" --run_only_null_model "
#null_model_str="  "

# default SNPID is chr:pos:ref_alt
# set option to output chr:pos:ref:alt
#snpid="  "
snpid=" --ref_colon_alt "

# default CHR name matches CHR in the VCF
# TOPMed imputation has chr21
# To remove "chr" to get 21 in the output add
#chrname=" --remove_chr_prefix "

# set min_maf
# GLI recommendation is 0.001
# Can set lower, say, 0.0005 as filtering also done at meta analysis
#
min_maf=" 0.001  "

# option to suppress output header
# a copy of the header : /data/GLI.header.txt
#output_header=" --suppress_output_header "
#output_header=" "

# output file names according to analysis plan
# 6 output files
# <file_prefix>.<trait>.<exposure>.<strata>.<file_suffix>.txt

# prefix to output file
file_prefix="PHASE2.INGI-FVG.EA.DBP"

# suffix to output file - upload date
# !! ADD CHR so don't clobbre
file_suffix="chr${chr}"

# !! robust standard errors were default in previous versions
rses=" --compute_robust_ses "
# parallel computing
num_mkl_threads=8
# --------------- DO NOT EDIT

covariate_str=" --covariates  ${covariates} "

dest=$(echo ${out_dir} | awk '{split($0,a,"/");print a[2]}')

if [[ ${dest} =~ "netapp" ]];then
echo "Analyses will be performed in a netapp storage. Therefore, output will not be compressed"
# commands
$mmap \
--trait \
${trait} \
${snpid} \
${rses} \
${exposure_str} \
${ExC_str} \
${exstrata_str} \
--g_by_exposure_analysis \
--min_maf \
0.001 \
${strata_str} \
${covariate_str} \
--file_prefix \
${file_prefix} \
--file_suffix \
${file_suffix} \
--phenotype_filename \
${pheno} \
--sparse_geno_analysis \
--tab_delimited \
--binary_genotype_filename \
${geno} \
--num_mkl_threads \
${num_mkl_threads} \
--phenotype_id \
CODPAZ \
--suppress_output_header
else
  echo "Analyses will not be performed in a netapp storage. Output will be compressed"

  $mmap \
  --trait \
  ${trait} \
  ${snpid} \
  ${rses} \
  ${exposure_str} \
  ${ExC_str} \
  ${exstrata_str} \
  --g_by_exposure_analysis \
  --min_maf \
  0.001 \
  ${strata_str} \
  ${covariate_str} \
  --file_prefix \
  ${file_prefix} \
  --file_suffix \
  ${file_suffix} \
  --phenotype_filename \
  ${pheno} \
  --sparse_geno_analysis \
  --tab_delimited \
  --binary_genotype_filename \
  ${geno} \
  --num_mkl_threads \
  ${num_mkl_threads} \
  --phenotype_id \
  CODPAZ \
  --suppress_output_header && gzip *.chr${chr}.txt # Compressing the analyses output
fi
