#!/usr/bin/env bash

#script to automate analysis of CHARGE's Depression-Blood Pressure study

model=$1
sex_group=$2

case $model in
  MODEL_one)
  #base=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one
  depr=$3
  BP=$4
  in_pheno=$5
  chr=$6
  fout=$7
  if [[ ${depr} == "qDEPR" && ${sex_group} == "M" ]]
  then
  cov="PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq AGExQDEP AGEsqxQDEP qDEPR"
  out_p=${fout}
  if [ -d "${out_p}" ]
  then
  echo "Output main directory already exists"
  else
    mkdir -p ${out_p} && echo "Output main directory doesn't exists. Creating Output main directory"
  fi
  if [[ ${BP} == "DP" ]]
  then
    compl_out=${out_p}
    if [ ! -d ${compl_out} ]
    then
    mkdir -p ${compl_out}
  else
    echo "Trait output folder already exists"
  fi
  logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/males/${depr}/${BP}
  if [ ! -d ${logs} ]
  then
  mkdir -p ${logs}
  else
    echo "Logs folder already exists"
  fi
  if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
  then
  echo "Output folder is not empty. Please check if analyses have already been done"
  exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
  echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
  echo "
  mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix males_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
  " | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
elif [[ ${BP} == "SP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out} ]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/males/${depr}/${BP}
if [ ! -d ${logs} ]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix males_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
elif [[ ${BP} == "PP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out} ]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/males/${depr}/${BP}
if [ ! -d ${logs} ]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait ${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix males_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=7G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
fi
  #FEMALES ONLY
elif [[ ${depr} == "qDEPR" && ${sex_group} == "F" ]]
then
cov="PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq AGExQDEP AGEsqxQDEP qDEPR"
out_p=${fout}
if [ -d "${out_p}" ]
then
echo "Output main directory already exists"
else
  mkdir -p ${out_p} && echo "Output main directory doesn't exists. Creating Output main directory"
fi
if [[ ${BP} == "DP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out}]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/female/${depr}/${BP}
if [ ! -d ${logs}]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N ${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
elif [[ ${BP} == "SP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out}]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/female/${depr}/${BP}
if [ ! -d ${logs}]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
elif [[ ${BP} == "PP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out}]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/female/${depr}/${BP}
if [ ! -d ${logs}]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait ${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
fi
elif [[ ${depr} == "qDEPR" && ${sex_group} == "MF" ]]
then
cov="PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq AGExQDEP AGEsqxQDEP MMAP_sex SEXxQDEP qDEPR"
out_p=${fout}
if [ -d "${out_p}" ]
then
echo "Output main directory already exists"
else
  mkdir -p ${out_p} && echo "Output main directory doesn't exists. Creating Output main directory"
fi
if [[ ${BP} == "DP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out}]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/sex_comb/${depr}/${BP}
if [ ! -d ${logs}]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
elif [[ ${BP} == "SP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out}]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/sex_comb/${depr}/${BP}
if [ ! -d ${logs}]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "f
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
elif [[ ${BP} == "PP" ]]
then
  compl_out=${out_p}
  if [ ! -d ${compl_out}]
  then
  mkdir -p ${compl_out}
else
  echo "Trait output folder already exists"
fi
logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_one/Logs/sex_comb/${depr}/${BP}
if [ ! -d ${logs}]
then
mkdir -p ${logs}
else
  echo "Logs folder already exists"
fi
if [[ ! -z "$(ls -A ${compl_out})" || $(ls -A ${compl_out}) > 0  ]]
then
echo "Output folder is not empty. Please check if analyses have already been done"
exit 1
elif [[ $(ls -A ${compl_out}) == 220 ]]
then
echo "Analyses for this BP trait are already completed"
exit 1
else
echo "Analyses files for ${depr} - ${BP} are not present. Starting analyses..."
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait ${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --gxe_interaction ${depr} --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}x${depr}_interaction_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP}_${depr} -j y -o /dev/null -V -l h_vmem=15G -wd ${compl_out} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
fi
fi

;;
MODEL_two)
#base=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_two
BP=$3
in_pheno=$4
chr=$5
fout=$6
logs=${fout}/Logs
#logs=/netapp05/analisi_nardone/CHARGE_DeprBP/MODEL_two/Logs
if [ ! -d ${fout} ]
then
mkdir -p ${fout} && echo "Base folder for model 2 doesn't exist. Creating a new one"
else
  echo "Base folder for Model 2 exists"
fi
if [ ! -d ${logs} ]
then
mkdir -p ${logs} && echo "Logs folder for model 2 doesn't exist. Creating a new one"
else
  echo "Logs folder for Model 2 exists"
fi

if [[ ${sex_group} == "M" ]]
then
out_p=${fout}
cov="PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq"
if [ ! -d ${out_p} ]
then
mkdir -p ${out_p} && echo "Folder for Model 2's analysis of ${BP} doesn't exist. Creating a new one"
else
  echo "Folder for Model 2's analysis of ${BP} already exists"
fi

spec_l=${logs}/males/${BP}

if [ ! -d ${spec_l} ]
then
mkdir -p ${spec_l} && echo "Creating logs folder for ${BP}"
else
  echo "Logs folder dor ${BP} already exists"
fi

if [[ ${BP} == "DP" || ${BP} == "SP" ]]
then
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix males_${BP}_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP} -j y -o /dev/null -V -l h_vmem=15G -wd ${out_p} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
elif [[ ${BP} == "PP" ]]
then
echo "f
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait ${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix males_${BP}_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP} -j y -o /dev/null -V -l h_vmem=15G -wd ${out_p} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi

elif [[ ${sex_group} == "F" ]]
then
out_p=${fout}/females_only/${BP}
cov="PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq"
if [ ! -d ${out_p} ]
then
mkdir -p ${out_p} && echo "Folder for Model 2's analysis of ${BP} doesn't exist. Creating a new one"
else
  echo "Folder for Model 2's analysis of ${BP} already exists"
fi

spec_l=${logs}/females/${BP}

if [ ! -d ${spec_l} ]
then
mkdir -p ${spec_l} && echo "Creating logs folder for ${BP}"
else
  echo "Logs folder dor ${BP} already exists"
fi

if [[ ${BP} == "DP" || ${BP} == "SP" ]]
then
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBPpheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP} -j y -o /dev/nullOB_ID_${BP}_${chr}.e -V -l h_vmem=15G -wd ${out_p} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
elif [[ ${BP} == "PP" ]]
then
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait ${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix females_${BP}_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP} -j y -o /dev/null -V -l h_vmem=15G -wd ${out_p} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi

elif [[ ${sex_group} == "MF" ]]
then
out_p=${fout}
cov="PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 AGE AGE_sq MMAP_sex"
if [ ! -d ${out_p} ]
then
mkdir -p ${out_p} && echo "Folder for Model 2's analysis of ${BP} doesn't exist. Creating a new one"
else
  echo "Folder for Model 2's analysis of ${BP} already exists"
fi

spec_l=${logs}/sex_comb/${BP}

if [ ! -d ${spec_l} ]
then
mkdir -p ${spec_l} && echo "Creating logs folder for ${BP}"
else
  echo "Logs folder dor ${BP} already exists"
fi

if [[ ${BP} == "DP" || ${BP} == "SP" ]]
then
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait Mean${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix comb_${BP}_${chr} --num_mkl_threads
" | qsub -N chr${chr}_${sex_group}_${BP} -j y -o /dev/null -V -l h_vmem=15G -wd ${out_p} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
elif [[ ${BP} == "PP" ]]
then
echo "
mmap --ped /netapp05/analisi_nardone/CHARGE_DeprBP/pheno_n_cov_files/FVG_Calcium_ECG_ped.csv --phenotype_filename ${in_pheno} --phenotype_id CODPAZ --trait ${BP} --covariate_filename ${in_pheno} --covariates ${cov} --read_binary_covariance_file /netapp04/concas/dose_MMAP/FVG/mmap_kinship.bin --single_pedigree --linear_regression --binary_genotype_filename /netapp04/concas/dose_MMAP/FVG/bin_files/chr${chr}_gen_binary --hc0_sandwich_estimator --file_suffix comb_${BP}_${chr} --num_mkl_threads 8
" | qsub -N chr${chr}_${sex_group}_${BP} -j y -o /dev/null -V -l h_vmem=15G -wd ${out_p} -q fast -pe smp 8 -m ea -M giuseppegiovanni.nardone@burlo.trieste.it
fi
fi
;;
esac
