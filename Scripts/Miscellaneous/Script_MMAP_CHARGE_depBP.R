#SEX combined dDEPR
#SPxdDEPR

for chr in {1..22}
do
echo "mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/FVG_Calcium_ECG_ped.csv --phenotype_filename /netapp05/analisi_nardone/CHARGE_DeprBP/FVG_Dep_BP_complete_sex_combined.csv --phenotype_id CODPAZ --trait MeanSP --covariate_filename /netapp05/analisi_nardone/CHARGE_DeprBP/FVG_Dep_BP_complete_sex_combined.csv --covariates PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq MMAP_sex AGExDDEP AGEsqxDDEP SEXxDDEP dDEPR --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction dDEPR --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix sex_comb_SPxdDEPR_interaction_${chr}" | qsub -N job${chr}_SPdDEP  -o /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/\$JOB_ID_${chr}_SPdDEP.log -e /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/\$JOB_ID_${chr}_SPdDEP.e -V -l h_vmem=7G -cwd -q fast -pe smp 8
done 

#DPxdDEPR

for chr in {1..22}
do
echo "mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Dep_BP_complete_sex_combined.csv --phenotype_id CODPAZ --trait MeanDP --covariate_filename /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Dep_BP_complete_sex_combined.csv --covariates PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq MMAP_sex AGExDDEP AGEsqxDDEP SEXxDDEP dDEPR --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction dDEPR --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix sex_comb_DPxdDEPR_interaction_${chr}" | qsub -N job${chr}_DPdDEP  -o /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/sex_comb/dDEPR/DP/\$JOB_ID_${chr}_DPdDEP.log -e /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/sex_comb/dDEPR/DP/\$JOB_ID_${chr}_DPdDEP.e -V -l h_vmem=7G -cwd -q fast -pe smp 8
done 

#PPxdDEPR
for chr in {1..22}
do
echo "mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Dep_BP_complete_sex_combined.csv --phenotype_id CODPAZ --trait PP --covariate_filename /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Dep_BP_complete_sex_combined.csv --covariates PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq MMAP_sex AGExDDEP AGEsqxDDEP SEXxDDEP dDEPR --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction dDEPR --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix sex_comb_DPxdDEPR_interaction_${chr}" | qsub -N job${chr}_DPdDEP  -o /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/sex_comb/dDEPR/PP/\$JOB_ID_${chr}_PPdDEP.log -e /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/sex_comb/dDEPR/PP/\$JOB_ID_${chr}_PPdDEP.e -V -l h_vmem=7G -cwd -q fast -pe smp 8
done 

#Males only
#SP
for chr in {1..22}
do
echo "mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_1927_complete_DepBP_males_only.csv --phenotype_id CODPAZ --trait MeanSP --covariate_filename /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_1927_complete_DepBP_males_only.csv --covariates PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq AGExDDEP AGEsqxDDEP dDEPR --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction dDEPR --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix sex_comb_DPxdDEPR_interaction_${chr}" | qsub -N job${chr}_SPdDEP  -o /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/males/SP/\$JOB_ID_${chr}_SPdDEP.log -e /netapp05/analisi_nardone/CHARGE_DeprBP/MODEL1/Logs/males/SP/\$JOB_ID_${chr}_SPdDEP.e -V -l h_vmem=7G -cwd -q fast -pe smp 8
done 
