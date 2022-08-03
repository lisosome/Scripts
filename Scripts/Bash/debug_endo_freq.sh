#!/usr/bin/env bash
# Author: Giuseppe Giovanni Nardone
#Aim: This script aims to locate variants in all imputed and genotyped populations. Moreover, extracts a list of carriers of the identified variations.
#The script was conceptualized for single variant research. However, it is possible to research regions. This version was ultimated on 08/02/2022.


#Display error message if one of the options is missing
if [[ $# -lt 4 ]]; then
  echo -e '\nError!!Missing arguments\n\n****** USAGE *****'
  echo -e 'Populations_variant_screening.sh insert options:'
  echo -e '-f <-- Specify a .bed input file.\n Bed file must contain chromosome, start and ending position as the 1st, 2nd and 3rd column respectively.\nFor single variants, start and ending positions must be the same!!'
  echo -e '-p <-- Choose a population. Options available: FVG, CAR, VBI, MOL [build 37], SR [Build 38]'
  echo -e '-b <-- Specify the build reference. Choose between GRCh37:"37",GRCh38: "38"'
  echo -e '-o <-- Specify an output folder'
  echo -e "All paths MUST BE ABSOLUTE!!\n"
   exit 1
fi
#getops to input the options, -f, -p, -b, -o.
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
#Core of the script. Depending onn the -b option, if the build condition is true, all above the else command will be executed
   if [[ ${build} == "37" ]]
     then
       p_name=(${paths37[${pop}]})
#echo "${p_name}"
       bcftools_bin=/share/apps/bio/bin/bcftools
       header_out="CHROM\tPOS\tID\tREF\tALT\tAC\tAN\tAF\tMAF\tAC_Het\tAC_Hom\tAC_Hemi\tCarriers"
       echo -e "${header_out}" > ${out}/Endometriosis_${pop}_frequency_test.txt
#While read scans the file line by line. For each file's line printed out by the echo, the three variables assume the values in the 1st, 2nd and 3rd column. In this case "lime" is the iterative variable of the loop. It can be named as you prefer (e.g Pippo, Pluto, Topolino).
         while read line
         do
           chr=$(echo -e $line | awk 'BEGIN{FS=" "} {print $1}')
           start=$(echo -e $line | awk 'BEGIN{FS=" "} {print $2}')
           end=$(echo -e $line | awk 'BEGIN{FS=" "} {print $3}')
#Principal command of the scripr. Bcftools view and query are used to locate regions and extract the desired info. Moreover, the plugin fill-tags is used for an easier extraction of information.
           (${bcftools_bin} view -r ${chr}:${start}-${end} ${p_name}/${chr}.vcf.gz | ${bcftools_bin} +fill-tags | ${bcftools_bin} query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\t%AF\t%MAF\t%AC_Het\t%AC_Hom\t%AC_Hemi[\t%SAMPLE=%GT]\n" | awk 'BEGIN{FS="\t"}{for(i=1;i<13;i++) {printf "%s ",$i };for(j=13;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",a[1]};printf "\n";}'| tr " " "\t")
           done < $inp_file >> ${out}/Endometriosis_${pop}_frequency_test.txt
    else
#If the starting conditionn is false, then this group of commands will be executed. They are exactly the same as above exept for for the bcftools view command: "chr" is added in order to locate regions in build 38.
       p_name=(${paths38[${pop}]})
       bcftools_bin=/share/apps/bio/bin/bcftools
       header_out="CHROM\tPOS\tID\tREF\tALT\tAC\tAN\tAF\tMAF\tAC_Het\tAC_Hom\tAC_Hemi\tCarriers"
       echo -e "${header_out}" > ${out}/Endometriosis_${pop}_frequency_test.txt
         while read line
        do
           chr=$(echo -e $line | awk 'BEGIN{FS=" "} {print $1}')
           start=$(echo -e $line | awk 'BEGIN{FS=" "} {print $2}')
           end=$(echo -e $line | awk 'BEGIN{FS=" "} {print $3}')
           (${bcftools_bin} view -r chr${chr}:${start}-${end} ${p_name}/${chr}.vcf.gz | ${bcftools_bin} +fill-tags | ${bcftools_bin} query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%AN\t%AF\t%MAF\t%AC_Het\t%AC_Hom\t%AC_Hemi[\t%SAMPLE=%GT]\n" | awk 'BEGIN{FS="\t"}{for(i=1;i<13;i++) {printf "%s ",$i };for(j=13;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",a[1]};printf "\n";}'| tr " " "\t")
           done < $inp_file >> ${out}/Endometriosis_${pop}_frequency_test.txt
       #echo 'Test passed'
  fi
