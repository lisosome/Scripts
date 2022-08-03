#!/usr/bin/env bash

#First we havo to generate the paths for our desired files. Provide a file with folder names and one with samples' names. Specify the date for traceability

smpl_folders=$1
infile=$2
trasf_path=$3
out_folder=$4
date=$5
field=$6

s3_path=s3://sftp.engenome.com/genmedburlo
readarray -t runs < ${smpl_folders}

for ((run_n=0;run_n<${#runs[@]};run_n++))
do
run_name=(${runs[${run_n}]})
  while read line
  do
    sam_n=$(echo -e $line)
    find -L /analisi_da_consegnare/burlo/${run_name} -type f -name "${sam_n}*.gz" | sort >> ${out_folder}/${date}_${run_name}_files_list.txt
  done < ${infile}
done

#The following loop is to set the correct names following illumina convention. Please verify that awk commands correctly split the file path
#out_path=/large/___SCRATCH___/burlo/nardone/220415_enG_transfer_files
#out_folder=/large/___HOME___/burlo/nardone/file_transfer

for ((n=0,x=1;n<${#runs[@]},x<=2;n++,x++))
do
  run_name=(${runs[${n}]})
  cat ${out_folder}/${date}_${run_name}_files_list.txt | awk -v outfile=${trasf_path} -v lane=${x} -v fld=${field} 'BEGIN{OFS=" "} {split($1,a,"/");split(a[fld],b,"_"); print outfile"/"b[1]"_"b[2]"_L00"lane"_"b[3]"_"b[4]}' >> ${out_folder}/${date}_transfer_lane${x}_namefixed_list.txt
done

#Generetion of executable files to perform the transfer in the scratch folder in order to create files with correctly renamed

for ((n=0,x=1;n<${#runs[@]},x<=2;n++,x++))
do
  run_name=(${runs[${n}]})
  (echo '#!/usr/bin/env bash';paste -d " " ${out_folder}/${date}_${run_name}_files_list.txt ${out_folder}/${date}_transfer_lane${x}_namefixed_list.txt | awk '{OFS=" "} {print "rsync","-avP",$0}') > ${out_folder}/${date}_lane_${x}_transfer.sh
done

#Transfer file to the $transf_path folder with fixed names. Then, execute md5

for n in 1 2
do
  chmod +x ${out_folder}/${date}_lane_${n}_transfer.sh
  source ${out_folder}/${date}_lane_${n}_transfer.sh
done

for ((run_n=0;run_n<${#runs[@]};run_n++))
do
run_name=(${runs[${run_n}]})
  while read line
  do
    sam_n=$(echo -e $line)
    find -L /analisi_da_consegnare/burlo/${run_name} -type f -name "${sam_n}*.gz" -exec md5sum {} \; | sort >> ${out_folder}/${date}_${run_name}_md5.txt
  done < ${infile}
done

for n in 1 2
do
find -L ${transf_path} -type f -name "*L00${n}*fastq.gz" -exec md5sum {} \; | sort >> ${out_folder}/${date}_lane${n}_renamed_md5.txt
done

for ((n=0,x=1;n<${#runs[@]},x<=2;n++,x++))
do
  run_name=(${runs[${n}]})
  diff <(cut -d " " -f 1 ${out_folder}/${date}_${run_name}_md5.txt | sort) <(cut -d " " -f 1 ${out_folder}/${date}_lane${x}_renamed_md5.txt | sort) > ${out_folder}/md5_check_lane${x}.txt
done

if [[ ! -s ${out_folder}/md5_check_lane1.txt && ! -s ${out_folder}/md5_check_lane2.txt ]]
then
  echo "Integrity check passed! Proceeding with the aws transfer"
  module load conda
  conda activate awscli
  aws s3 cp ${trasf_path}/ ${s3_path}/${date}/ --exclude "*" --include "*.fastq.gz"  --recursive --acl bucket-owner-full-control
  exit 0
else
  echo " Data integrity corrupted. Cleaning All"
  rm ${transf_path}/*
  exit 1
fi
