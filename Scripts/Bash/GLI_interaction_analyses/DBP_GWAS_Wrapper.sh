#!/usr/bin/env bash

# Wrapper script for easier job submission for CHARGE Depression-BP GWAS analyses

# Setting variables. They will be the same as the GWAS script, with the add of the option for the output folder
mode=$1

case $mode in
  SINGLE_CHR)
  chr=$2
  trt=$3
  out_dir=$4

echo "SINGLE CHROMOSOME MODE:
       chromosome: $chr
       trait: $trt
       output directory: $out_dir"

if [[ ! -d ${out_dir} ]]; then
  echo -e "\nMissing output directory. Creating the output directory and a log subfolder"
  mkdir -p ${out_dir}/Logs ${out_dir}/csv_files
  echo "$out_dir, $out_dir/Logs and ${out_dir}/csv_files created. Proceeding with the analyses for th $trt trait"
elif [[ ! -d ${out_dir}/Logs && -d ${out_dir}/csv_files ]];then
echo -e "\nMissing Log folder! Creating a new one"
mkdir -p ${out_dir}/Logs
elif [[ -d ${out_dir}/Logs && ! -d ${out_dir}/csv_files ]];then
echo -e "\nMissing folder for csv files! Creating a new one"
mkdir -p ${out_dir}/csv_files
else
  echo -e "\nExisting output directory. \nProceeding with the analyses for the $trt trait "
fi

qsub -N BDP_GWAS_${chr} -q fast -pe smp 8 -wd ${out_dir} -o ${out_dir}/Logs/DBP.GWAS_chr${chr}.log -e ${out_dir}/Logs/DBP.GWAS_chr${chr}.e /home/nardone/script/CHARGE_miniscript/DBP_GWAS.sh $chr $trt $out_dir &&
qsub -N chr${chr}_check -q fast -pe smp 8 -wd ${out_dir} -j y -o /dev/null -hold_jid "BDP_GWAS_${chr}" /home/nardone/script/CHARGE_miniscript/DBP_analyses_check.sh $chr $trt $out_dir

;;

  ALL_CHR)
  trt=$2
  out_dir=$3

  echo "ALL_CHR mode selected. Jobs for autosomal chromosomes (1-22) will be submitted for $trt trait
        trait: $trt
        output directory: $out_dir"

  if [[ ! -d ${out_dir} ]]; then
    echo "Missing output directory. Creating the output directory and a log subfolder"
    mkdir -p ${out_dir}/Logs ${out_dir}/csv_files
    echo "$out_dir, $out_dir/Logs and ${out_dir}/csv_files created. Proceeding with the analyses for th $trt trait"
  elif [[ ! -d ${out_dir}/Logs && -d ${out_dir}/csv_files ]];then
  echo "Missing Log folder! Creating a new one"
  mkdir -p ${out_dir}/Logs
  elif [[ -d ${out_dir}/Logs && ! -d ${out_dir}/csv_files ]];then
  echo "Missing folder for csv files! Creating a new one"
  mkdir -p ${out_dir}/csv_files
  else
    echo "Existing output directory. Proceeding with the analyses for the $trt trait "
  fi

for cro in {1..22};do
qsub -N BDP_GWAS_${cro} -q fast -pe smp 8 -wd ${out_dir} -o ${out_dir}/Logs/DBP.GWAS_chr${cro}.log -e ${out_dir}/Logs/DBP.GWAS_chr${cro}.e /home/nardone/script/CHARGE_miniscript/DBP_GWAS.sh $cro $trt $out_dir &&
qsub -N chr$cro_check -q fast -pe smp 8 -wd ${out_dir} -j y -o /dev/null -hold_jid "BDP_GWAS_${cro}" /home/nardone/script/CHARGE_miniscript/DBP_analyses_check.sh $cro $trt $out_dir
done

;;
esac
