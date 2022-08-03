#!/usr/bin/env bash

panel_mode=$1
in_bed=$2
out_name=$3
build=$4
rec=$5

bedt=/engenome/codebase/WGA/tools/BEDTools-2.22.0/bedtools
out_folder=/engenome/libraries/${build}/${out_name}_panel
anno37=/engenome/codebase/WGA/bundle/current/Ensembl75_annotation_track_v3/Ensemble75_v37.out.sorted.bed
anno38=/engenome/codebase/WGA/bundle/hg38/Ensembl97_annotation_track_v2/HG38/Ensemble97_v38.HG38.out.sorted.bed

if [[ $build == "GRCh37" || $build == "hg38" ]]
then
  mkdir -p ${out_folder}
else
  echo "Invalid genome reference. Please choose between GRCh37 and hg38"
exit 1
fi

echo "Library folder ${out_folder} created"

if [[ $panel_mode == "custom" ]]
then
  if [[ $build == "GRCh37" ]]
  then
    awk -v rd=${rec} '(NR >= rd){print substr($0,4)}' ${in_bed} | ${bedt} sort -i - | ${bedt} merge -i - | ${bedt} sort -i - > ${out_folder}/${out_name}_panel.bed &&
    echo "Cleaned library bed file available at ${out_folder}/${out_name}"
elif [[  $build == "hg38" ]]
then
  awk -v rd=${rec} '(NR >= rd)' ${in_bed} | ${bedt} sort -i - | ${bedt} merge -i - | ${bedt} sort -i - > ${out_folder}/${out_name}_panel.bed &&
  echo "Cleaned library bed file available at ${out_folder}/${out_name}"
else
  echo "Please select a valid genome reference. Choose between GRCh37 and hg38"
exit 1
fi
elif [[ $panel_mode == "commercial" ]]
then
  if [[ $build == "GRCh37" ]]
  then
    cat $in_bed | grep -v rand | grep -v Un | grep -v M | grep -v hap | bedtools sort -i - | bedtools merge -i - | grep -v “^\#” | bedtools sort -i - > ${out_folder}/${out_name}_panel.bed &&
    echo "Cleaned library bed file available at ${out_folder}/${out_name}"
elif [[ $build == "hg38" ]]
then
  cat $in_bed | grep -v rand | grep -v Un | grep -v M | grep -v hap | bedtools sort -i - | bedtools merge -i - | grep -v “^\#” | bedtools sort -i - > ${out_folder}/${out_name}_panel.bed &&
  echo "Cleaned library bed file available at ${out_folder}/${out_name}"
fi
fi

### Gene list generation ###

if [[ $build == "GRCh37" ]]
then
  bash /engenome/codebase/WGA/tools/pipeline_germline_engenome_v2.3/scripts/Build_gene_list_v0.1.sh -b ${out_folder}/${out_name}_panel.bed -g ${anno37} -p ${out_folder} -o gene_list -l ${bedt}
### CNV bed generation ###
mkdir -p ${out_folder}/CNV/controls/Females ${out_folder}/CNV/controls/Males
bash /engenome/codebase/WGA/tools/pipeline_germline_engenome_v2.3/scripts/Prepare_bedv4.3.sh -m 1 -t ensembl -b ${out_folder}/${out_name}_panel.bed -f 2 -g ${anno37} -p ${out_folder}/CNV -o ${out_name} -l ${bedt} -s /engenome/codebase/WGA/tools/singularity_v3.0.3/workspace/go/bin/bin/singularity -u /engenome/codebase/WGA/tools/singularity_v3.0.3/containers/ubuntu_18.04.sif
mv ${out_folder}/CNV/${out_name}.complete.bed ${out_folder}/CNV/${out_name}.CNV.bed
elif [[ $build == "hg38" ]]
then
  bash /engenome/codebase/WGA/tools/pipeline_germline_engenome_v2.3/scripts/Build_gene_list_v0.1.sh -b ${out_folder}/${out_name}_panel.bed -g ${anno38} -p ${out_folder} -o gene_list -l ${bedt}
### CNV bed generation ###
mkdir -p ${out_folder}/CNV/controls/Females ${out_folder}/CNV/controls/Males
bash /engenome/codebase/WGA/tools/pipeline_germline_engenome_v2.3/scripts/Prepare_bedv4.3.sh -m 1 -t ensembl -b ${out_folder}/${out_name}_panel.bed -f 2 -g ${anno38} -p ${out_folder}/CNV -o ${out_name} -l ${bedt} -s /engenome/codebase/WGA/tools/singularity_v3.0.3/workspace/go/bin/bin/singularity -u /engenome/codebase/WGA/tools/singularity_v3.0.3/containers/ubuntu_18.04.sif
mv ${out_folder}/CNV/${out_name}.complete.bed ${out_folder}/CNV/${out_name}.CNV.bed
fi
