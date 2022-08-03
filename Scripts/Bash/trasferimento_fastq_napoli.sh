#!/usr/bin/env bash

dest_path=/net/backup/TRANSFER_DISTRO_RET
inp_file=$1
new_fol=$2

finfold=${dest_path}/${new_fol}

while read line
do
  base_path=$(echo -e $line | awk 'BEGIN{FS=" "} print $2}')
  sample=$(echo -e $line | awk 'BEGIN{FS=" "} print $1}')
  #echo "$base_path"
  #echo "$sample"
  mkdir -p ${finfold}/${sample}
  rsync -avP ${base_path}/${sample}*.gz ${finfold}/${sample}/.
  find ${base_path} -type f -name "${sample}*.gz" -exec md5sum {} \; >> ${HOME}/md5_original
  find ${finfold}/${sample} -type f -name "${sample}*.gz" -exec md5sum {} \; >> ${HOME}/md5_transferred
done < ${inp_file}

check_transf=$(diff <(cut -d " " -f1 ${HOME}/md5_original | sort ) <(cut -d " " -f1 ${HOME}/md5_original | sort))

if [[ -z ${check_transf} ]]
then
  echo "Transfer completed succesfully"
  rm ${HOME}/md5_original ${HOME}/md5_transferred
else
  echo "Transfer failed. Something went wrong"
  rm -r ${finfold}
fi
