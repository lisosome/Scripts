#!/usr/bin/env bash

declare -a LOF=(splice_acceptor_variant splice_donor_variant frameshift_variant stop_lost start_lost stop_gained missense_variant NMD_transcript_variant)
header="CHROM\tPOS\tID\tREF\tALT\tAF\tAC_Hemi\tGene\tConsequence\tCADD_PHRED\tTranscript/Total_transcripts\tTranscript_Score\tCarriers"
base_path=/home/nardone/COVID/Lof/retake/retake2
mkdir -p ${base_path}

for lof in ${!LOF[@]}
do
  var=${LOF[${lof}]}
  file_out=${base_path}/${var}_AF_noZero_CADD.txt
  echo -e "${header}" > ${file_out}
  for chr in {1..22} X
  do
    bcftools view -e "INFO/AF=='0'" /netapp06/WGS/RELEASE/20220210/GRCh38/OTORINO/1.PHASED/${chr}.vcf.gz | /share/apps/bio/bin/bcftools +split-vep -i "Consequence~'${var}'" -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AC_Hemi\t%SYMBOL\t%Consequence\t%CADD_PHRED[\t%SAMPLE=%GT]\n" | awk 'BEGIN{FS="\t";OFS="\t"} {split($8,x,","); $8=x[1]} {print $0}' |awk 'BEGIN{FS="\t";OFS="\t"} {split($10,x,","); $10=x[1]} {print $0}'| awk -v lof=$var '{FS="\t";OFS="\t"}{cont=0;n=split($9,a,",");{for(i in a);if(a[i]=="lof");cont=cont+1};$11=cont"/"n;{print $0}}' | awk -v lof=$var '{FS="\t";OFS="\t"}{cont=0;n=split($9,a,",");{for(i in a);if(a[i]=="lof");cont=cont+1};$12=cont/n;{print $0}}' | awk 'BEGIN{FS="\t"}{for(i=1;i<13;i++) {printf "%s ",$i };for(j=13;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",$j};printf "\n";}'| tr " " "\t"  >> ${file_out}
  done
done




#for lof in ${LOF[@]}
#do
  #file_out=${base_path}/${lof}_AF_noZero_CADD.txt
  #cat ${file_out} | awk -v lof=$lof '{FS="\t";OFS="\t"}
    #{cont=0
    #n=split($8,a,",")
    #{for (i in a)
      #if(a[i]==lof)
    #cont=cont+1}
    #$10=print cont"/"n
    #{print $0}
  #}'

  #awk -v lof=$lof '{FS="\t";OFS="\t"}{cont=0;n=split($8,a,",");{for (i in a);if(a[i]==lof) cont=cont+1};$10=print cont"/"n;{print $0}}'
