#!/usr/bin/env bash
# Author: Giuseppe Giovanni Nardone
#Aim: This script aims to locate variants in all imputed and genotyped populations. Moreover, extracts a list of carriers of the identified variations.
#The script was conceptualized for single variant research. However, it is possible to research regions. This version was ultimated on 08/02/2022.


#Display error message if one of the options is missing
if [[ $# -lt 1 ]]; then
  echo -e '\nError!!Missing arguments\n\n****** USAGE *****'
  echo -e 'Populations_variant_screening.sh insert options:\n'
  echo -e '-f <-- Specify a .bed input file.\n\tBed file must contain chromosome (e.g. 1 for 37 or chr1 for 38), start and ending position as the 1st, 2nd and 3rd column respectively.\n\tFor single variants, start and ending positions must be the same!!'
  echo -e '-p <-- Choose a population. Options available: FVG, CAR, VBI, MOL [build 37], SR, HCFVG [Build 38]'
  echo -e '-b <-- Specify the build reference. Choose between GRCh37:"37",GRCh38: "38"'
  echo -e '-o <-- Specify an output folder\n'
  echo -e "\tAll paths MUST BE ABSOLUTE!!\n"
   exit 1
fi

#Display help message
if [[ $@ == "-h" ]]
then
echo -e '\nPopulations_variant_screening.sh insert options:\n'
echo -e '-f <-- Specify a .bed input file.\n\tBed file must contain chromosome (e.g. 1 for 37 or chr1 for 38), start and ending position as the 1st, 2nd and 3rd column respectively.\n\tFor single variants, start and ending positions must be the same!!'
echo -e '-p <-- Choose a population. Options available: FVG, CAR, VBI, MOL [build 37], SR, HCFVG [Build 38]'
echo -e '-b <-- Specify the build reference. Choose between GRCh37:"37",GRCh38: "38"'
echo -e '-o <-- Specify an output folder\n'
echo -e "\tAll paths MUST BE ABSOLUTE!!\n"
exit 1
fi

#getops to input the options, -f, -p, -b, -o,
echo "${@}"
while getopts ":f:p:b:o:" opt ${@}; do
        case $opt in
          f)
             echo ${OPTARG}
             inp_file=${OPTARG}
          ;;
          p)
             echo ${OPTARG}
             pop=${OPTARG}
          ;;
          b)
             echo ${OPTARG}
             build=${OPTARG}
          ;;
          o)
             echo ${OPTARG}
             out=${OPTARG}
          ;;
       esac
done
#Creation of two associative arrays collecting paths for genotyped and imputed data
declare -A paths37=()
paths37[FVG]=$(echo '/shared/INGI_WGS/20200708/FVG')
paths37[CAR]=$(echo '/shared/INGI_WGS/20200708/CAR')
paths37[VBI]=$(echo '/shared/INGI_WGS/20200708/VBI')
paths37[MOL]=$(echo "/netapp02/data/imputation/IGRPv1/MOLISANI/20210613/06.imputed/VCF/${chr}")

declare -A paths38=()
paths38[SR]=$(echo '/shared/WGS/SilkRoad/20210322_RELEASE_FILTERED')
paths38[HCFVG]=$(echo '/netapp06/WGS/RELEASE/20220210/GRCh38/FVG/1.PHASED')

#Check that every option is set correctly
if [[ ${build} == "38" && (${pop} == "FVG" || ${pop} == "CAR" || ${pop} == "VBI" || ${pop} == "MOL") ]]
then
  echo -e '***ERROR***\n\nPlease check the correct build version for the selected population. FVG, CAR, VBI, MOL [build 37], SR, HCFVG [Build 38]\n'
exit 1
elif [[ ${build} == "37" && (${pop} == "SR" || ${pop} == "HC_FVG") ]]
then
  echo -e '***ERROR***\n\nPlease check the correct build version for the selected population. FVG, CAR, VBI, MOL [build 37], SR, HCFVG [Build 38]\n'
exit 1
fi


#Core of the script. Depending on the -b option, if the build condition is true, all above the else command will be executed
if [[ ${build} == "37" ]]
     then
       p_name=(${paths37[${pop}]})
#echo "${p_name}"
       bcftools_bin=/share/apps/bio/bin/bcftools
       header_out="CHROM\tPOS\tID\tREF\tALT\tAC\tAN\tAF\tMAF\tAC_Het\tAC_Hom\tAC_Hemi\tGene\tIsoform\tConsequence\tgnomAD_AF\tCDS_position\tExon\tIntron\tProtein_position\tAmino_acids\tSIFT\tPolyPhen\tCADD_RAW\tCarriers"
       echo -e "${header_out}" > ${out}/${pop}.txt
#While read scans the file line by line. For each file's line printed out by the echo, the variable assume the values in the 1st column. In this case "line" is the iterative variable of the loop. It can be named as you prefer (e.g Pippo, Pluto, Topolino).
         while read line
         do
          chr=$(echo -e $line | awk 'BEGIN{FS=" "} {sub(/chr/, ""); print $1}')
          ${bcftools_bin} view -R ${inp_file} ${p_name}/${chr}.vcf.gz | ${bcftools_bin} +fill-tags | ${bcftools_bin} query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\t%AF\t%MAF\t%AC_Het\t%AC_Hom\t%AC_Hemi\n" >> fill_tags.txt
          ${bcftools_bin} view -R ${inp_file} ${p_name}/${chr}.vcf.gz | ${bcftools_bin} +split-vep -d -f "%ID\t%SYMBOL\t%Feature\t%Consequence\t%gnomAD_AF\t%CDS_position\t%EXON\t%INTRON\t%Protein_position\t%Amino_acids\t%SIFT\t%PolyPhen\t%CADD_RAW[\t%SAMPLE=%GT]\n" | awk 'BEGIN{FS="\t"}{for(i=1;i<14;i++) {printf "%s ",$i };for(j=14;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",a[1]};printf "\n";}'| tr " " "\t" >> split_vep.txt
         done < ${inp_file}
        join -1 3 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 fill_tags.txt split_vep.txt > ${out}/${pop}_variants_screening.txt
        awk '!seen[$2,$3,$15]++' ${out}/${pop}_variants_screening.txt  >> ${out}/${pop}.txt
        rm fill_tags.txt split_vep.txt ${out}/${pop}_variants_screening.txt
else
#If the starting condition is false, then this group of commands will be executed. They are exactly the same as above exept for for the chr variable.
       p_name=(${paths38[${pop}]})
       bcftools_bin=/share/apps/bio/bin/bcftools
       header_out="CHROM\tPOS\tID\tREF\tALT\tAC\tAN\tAF\tMAF\tAC_Het\tAC_Hom\tAC_Hemi\tGene\tIsoform\tConsequence\tgnomAD_AF\tCDS_position\tExon\tIntron\tProtein_position\tAmino_acids\tSIFT\tPolyPhen\tCADD_RAW\tCarriers"
       echo -e "${header_out}" > ${out}/${pop}.txt
       while read line
       do
        chr=$(echo -e $line | awk 'BEGIN{FS=" "} {sub(/chr/, ""); print $1}')
        ${bcftools_bin} view -R ${inp_file} ${p_name}/${chr}.vcf.gz | ${bcftools_bin} +fill-tags | ${bcftools_bin} query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\t%AF\t%MAF\t%AC_Het\t%AC_Hom\t%AC_Hemi\n" >> fill_tags.txt
        ${bcftools_bin} view -R ${inp_file} ${p_name}/${chr}.vcf.gz | ${bcftools_bin} +split-vep -d -f "%CHROM\t%POS\t%ID\t%SYMBOL\t%Feature\t%Consequence\t%gnomAD_AF\t%CDS_position\t%EXON\t%INTRON\t%Protein_position\t%Amino_acids\t%SIFT\t%PolyPhen\t%CADD_RAW[\t%SAMPLE=%GT]\n" | awk 'BEGIN{FS="\t"}{for(i=1;i<16;i++) {printf "%s ",$i };for(j=16;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",a[1]};printf "\n";}'| tr " " "\t" >> split_vep.txt
       done < ${inp_file}
       #sort -k2,2n split_vep.txt > sorted_split.txt
       #sort -k2,2n fill_tags.txt > sorted_fill.txt
      join -1 3 -2 3 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14,2.15,2.16 fill_tags.txt split_vep.txt > ${out}/${pop}_variants_screening.txt
      awk '!seen[$2,$3,$15]++' ${out}/${pop}_variants_screening.txt  >> ${out}/${pop}.txt
      rm fill_tags.txt split_vep.txt ${out}/${pop}_variants_screening.txt
fi
