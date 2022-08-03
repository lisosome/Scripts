#!/usr/bin/env bash

# Script to concatenate all GLI interaction results file

base=$1
trait=$2
exp=$3
sex=$4

head=/home/nardone/software/mmap/mmap_gli_training/data/GLI.header.txt

if [[ $trait == "DP" ]];then filename=DBP
elif [[ $trait == "SP" ]];then filename=SBP
elif [[ $trait == "PPr" ]];then filename=PP
else
  echo "Select a valid trait. Choose between DP SP or PP"
fi

cat ${head} ${base}/${trait}/${exp}/${sex}/PHASE2.INGI-FVG.EA.DBP.${trait}.${exp}.${sex}.chr{1..22}.txt | sed 's/-nan/./g' |sed 's/nan/./g' | sed 's/inf/./g' >> ${base}/Final_files/${exp}/PHASE2.INGI-FVG.EA.${filename}.${exp}.${sex}.20220802.txt &&
gzip ${base}/Final_files/${exp}/PHASE2.INGI-FVG.EA.${filename}.${exp}.${sex}.20220802.txt
