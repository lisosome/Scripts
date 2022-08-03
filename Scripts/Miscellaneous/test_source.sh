#!/usr/bin/env bash
for chr in {1..22} X Y
do
source /large/___HOME___/burlo/nardone/script/WGS_various_analyses_106.sh CSQANNOT /large/___HOME___/burlo/nardone/hearing_cnv/surv_merge_sorted.vcf /large/___HOME___/burlo/nardone/hearing_cnv/chr${chr}_tst_merged_annotated_.vcf.gz ${chr} NEW ORFEO 106 GRCh38 /fast/burlo/fcrudele/resources/hgRef/hg38/Homo_sapiens_assembly38.fasta
done
