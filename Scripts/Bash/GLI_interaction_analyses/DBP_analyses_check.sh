#!/usr/bin/env bash

# Script to check if all files were generated during the analyses run

chr=$1
trt=$2
out_dir=$3

  if [[ -s ${out_dir}/Logs/DBP.GWAS_chr${chr}.e || ! -f ${out_dir}/*.chr${chr}*.gz  ]]; then
    echo "$trt trait association results for chromosome $chr are missing/empty!!
          Re-excecuting analyses for chromosome $chr
          Re-execution in 60 seconds..."
          sleep 1m
          qsub -N re_BDP_GWAS_${chr} -q fast -pe smp 8 -wd ${out_dir} -o ${out_dir}/Logs/DBP.GWAS_chr${chr}.log -e ${out_dir}/Logs/DBP.GWAS_chr${chr}.e /home/nardone/script/CHARGE_miniscript/DBP_GWAS.sh $chr $trt $out_dir
  else
    dest=$(realpath ${out_dir}/*.chr${chr}*.gz)
    echo "Checking for analyses files for chr$chr...
          $trt trait association results for chromosome $chr are available at:
          ${dest}"
  fi
