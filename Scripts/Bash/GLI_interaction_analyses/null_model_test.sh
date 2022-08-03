#!/usr/bin/env bash 

# Step 1 for related s script to run null model in MMAP
tra=$1
#mmap=../programs/mmap

# pedigree file  
ped=/netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv

# phenotype file
pheno=/netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Dep_BP_complete_sex_combined.csv

# !! SET CHOICE OF COVARIANCE MATRIX
# kinship from pedigree
covariance_matrix=/netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin

# grm from genotypes
#covariance_matrix=../grm2mmap/plink/5K.500K.simulated.rel.zstd.bin
#covariance_matrix=../grm2mmap/plink/5K.500K.simulated.rel.bin
#covariance_matrix=../grm2mmap/subject_included/5K.subject.included.bin
#covariance_matrix=../grm2mmap/subject_included/5K.subject.included.gz.bin

# trait and covariates assumed to be in phenotype file
trait=${tra}

# num threads, increase for faster compute/large samples
num_mkl_threads=8

# CHOOSE COVARIATES TO ADJUST OUT
covariates=" SEX AGE AGE_sq   " 
covariate_str=" --covariates ${covariates} " 

# if NO covariate adjustment, remove comment 
#covariate_str="  " 

# set suffix for MMAP null model files 
#mmap_output_suff=kinship
mmap_output_suff=null.mod


mmap \
--ped \
${ped} \
--trait \
$trait \
${covariate_str} \
--file_suffix \
$mmap_output_suff \
--phenotype_filename \
$pheno \
--single_pedigree \
--num_mkl_threads ${num_mkl_threads} \
--read_binary_covariance_file \
${covariance_matrix} 


