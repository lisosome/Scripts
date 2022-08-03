#set -e
#mode=$1
#case $mode in
#  VCFAnn)
if [[ $@ == "-h" ]]
then
  echo -e 'Annotate your VCF files with VEP!\n\nPlease insert the following options:'
  echo -e '-i <-- Input folder of your vcf files. The files must be named "chr.vcf.gz" (e.g. 1.vcf.gz, 2.vcf.gz)'
  echo -e '-b <-- Select your build of reference (GRCh37 or GRCh38)'
  echo -e '-o <-- Output folder of your annotated vcf files'
  echo -e '-v <-- Please select your preferred VEP software version (105, 100)'
  exit 1
fi
#To implement:other version of vep (90,94,105)
echo "${@}"
while getopts ":i:b:o:v:" opt ${@}; do
        case $opt in
          i)
             echo ${OPTARG}
             infolder=${OPTARG}
          ;;
          b)
             echo ${OPTARG}
             build=${OPTARG}
          ;;
          o)
             echo ${OPTARG}
             outfolder=${OPTARG}
          ;;
          v)
             echo ${OPTARG}
             vep_version=${OPTARG}
          ;;
       esac
done

maincache=/shared/resources/VEP_cache
gerp_db=${maincache}/${vep_version}/compara/gerp_conservation_scores.homo_sapiens.${build}.bw
dbNSFP_fields="Ensembl_transcriptid,Uniprot_acc,VEP_canonical,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,
Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,
MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_
rankscore,MetaLR_pred,REVEL_score,REVEL_rankscore,ClinPred_score,ClinPred_rankscore,ClinPred_pred"
gnomad_annot_string="AC,AN,AF,nhomalt,AC-XY,AN-XY,AF-XY,nhomalt-XY,AC-oth,AN-oth,AF-oth,nhomalt-oth,AC-ami,AN-ami,AF-ami,nhomalt-ami,AC-sas,AN-sas,AF-sas,nhomal
t-sas,AC-fin,AN-fin,AF-fin,nhomalt-fin,AC-eas,AN-eas,AF-eas,nhomalt-eas,AC-amr,AN-amr,AF-amr,nhomalt-amr,AC-afr,AN-afr,AF-afr,nhomalt-afr,AC-mid,AN-mid,AF-mid,n
homalt-mid,AC-asj,AN-asj,AF-asj,nhomalt-asj,AC-nfe,AN-nfe,AF-nfe,nhomalt-nfe"
declare -a CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

case ${build} in
  GRCh38)
  cadd_v="v1.6"
  gnomad_v=${maincache}/${vep_version}/cadd/${cadd_v}/${build}/gnomad.genomes.r3.0.snv.tsv.gz
  gnomad_indel=${maincache}/${vep_version}/cadd/${cadd_v}/${build}/gnomad.genomes.r3.0.indel.tsv.gz
  cadd_whole_genome=${maincache}/${vep_version}/cadd/${cadd_v}/${build}/whole_genome_SNVs.tsv.gz
if [[ ${vep_version} == "105" ]]
then
  source activate vep105
  plugin_path=/home/santin/.conda/envs/vep105/share/ensembl-vep-105.0-0
  vep_bin=/home/santin/.conda/envs/vep105/share/ensembl-vep-105.0-0/vep
  for chr in ${CHR[@]}
  do
    gnomad_vcf=/netapp05/analisi_nardone/gnomAD/GRCh38/gnomad.genomes.r3.0.sites.chr${chr}_trimmed_info.vcf.bgz
    dbNSFP=/netapp05/analisi_nardone/dbNSFP/${build}/dbNSFP4.3a_grch38_chr${chr}.gz
    ${vep_bin} -i ${infolder}/${chr}.vcf.gz --stats_text --force_overwrite --format vcf --offline --fasta /shared/resources/hgRef/hg38/Homo_sapiens_assembly38.fasta --compress_output bgzip --dir_plugins ${plugin_path} --o ${outfolder}/vep.${chr}.vcf.gz --buffer_size 200000 --fork 12 --plugin CADD,${cadd_whole_genome},${gnomad_indel},${gnomad_v} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${build} --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin Conservation,${gerp_db} --plugin ExACpLI,${plugin_path}/ExACpLI_values.txt
done
  else
    source activate vep
    plugin_path=/home/nardone/.conda/envs/vep/share/ensembl-vep-100.2-0
    vep_bin=/home/nardone/.conda/envs/vep/share/ensembl-vep-100.2-0/vep
    for chr in ${CHR[@]}
    do
      gnomad_vcf=/netapp05/analisi_nardone/gnomAD/GRCh38/gnomad.genomes.r3.0.sites.chr${chr}_trimmed_info.vcf.bgz
      dbNSFP=/netapp05/analisi_nardone/dbNSFP/${build}/dbNSFP_4.1_variant.chr${chr}.gz
      ${vep_bin} -i ${infolder}/${chr}.vcf.gz --stats_text --force_overwrite --format vcf --offline --fasta /shared/resources/hgRef/hg38/Homo_sapiens_assembly38
.fasta --compress_output bgzip --dir_plugins ${plugin_path} --o ${outfolder}/vep.${chr}.vcf.gz --buffer_size 200000 --plugin CADD,${cadd_whole_genome},${gnomad_
indel},${gnomad_v} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds -
-protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --asse
mbly ${build} --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_posi
tion:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${plugin_path} --plugin Conservation,${gerp
_db}
  done
fi
;;
  GRCh37)
  cadd_v="v1.6"
  gnomad_v=${maincache}/${vep_version}/cadd/${cadd_v}/${build}/gnomad.genomes.r2.1.1.snv.tsv.gz
  gnomad_indel=${maincache}/${vep_version}/cadd/${cadd_v}/${build}/gnomad.genomes.r2.1.1.indel.tsv.gz
  cadd_whole_genome=${maincache}/${vep_version}/cadd/${cadd_v}/${build}/whole_genome_SNVs.tsv.gz
  if [[ ${vep_version} == "105" ]]
  then
    source activate vep105
    plugin_path=/home/santin/.conda/envs/vep105/share/ensembl-vep-105.0-0
    vep_bin=/home/santin/.conda/envs/vep105/share/ensembl-vep-105.0-0/vep
  for chr in ${CHR[@]}
  do
gnomad_vcf=/netapp05/analisi_nardone/gnomAD/GRCh37/gnomad.genomes.r2.1.sites.chr${chr}_noVEP.vcf.gz
dbNSFP=/netapp05/analisi_nardone/dbNSFP/GRCh37/dbNSFP4.3a_grch37_${chr}.gz
${vep_bin} -i ${infolder}/${chr}.vcf.gz --stats_text --force_overwrite --format vcf --offline --fasta /shared/resources/gatk4hg19db/Homo_sapiens_assembly19_1000genomes_decoy.fasta --compress_output bgzip --dir_plugins ${plugin_path} --o ${outfolder}/vep.${chr}.vcf.gz --buffer_size 200000 --fork 12 --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin CADD,${cadd_whole_genome},${gnomad_indel},${gnomad_v} --plugin Conservation,${gerp_db} --plugin ExACpLI,${plugin_path}/ExACpLI_values.txt --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${build}
done
 else
  source activate vep
  plugin_path=/home/santin/.conda/envs/vep/share/ensembl-vep-100.2-0
  vep_bin=/home/santin/.conda/envs/vep/share/ensembl-vep-100.2-0/vep
  for chr in ${CHR[@]}
  do
    gnomad_vcf=/netapp05/analisi_nardone/gnomAD/GRCh37/gnomad.genomes.r2.1.sites.chr${chr}_noVEP.vcf.gz
    dbNSFP=/netapp05/analisi_nardone/dbNSFP/dbNSFP_4.1_variant.chr${chr}.gz
    ${vep_bin} -i ${infolder}/${chr}.vcf.gz --stats_text --force_overwrite --format vcf --offline --fasta Homo_sapiens_assembly19_1000genomes_decoy.fasta --comp
ress_output bgzip --dir_plugins ${plugin_path} --o ${outfolder}/vep.${chr}.vcf.gz --buffer_size 200000 --plugin CADD,${cadd_whole_genome},${gnomad_indel},${gnom
ad_v} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --un
iprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${build}
 --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min
_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${plugin_path} --plugin Conservation,${gerp_db}
done
fi
;;
esac
#--plugin CADD,${cadd_whole_genome},${gnomad_indel},${gnomad_v} --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin LoF,human_ancestor_fa:${maincache}/${vep_ver
sion}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,lofte
e_path:${plugin_path} --plugin Conservation,${gerp_db}
#;;
#GWAS_BEST)
#if [[ $@ == "-h" ]]
#then
#  echo -e 'Annotate your Best SNPs file with VEP!\n\nPlease insert the following options:'
#  echo -e '-i <-- Absolute path of the input file'
#  echo -e '-b <-- Select your build of reference (GRCh37 or GRCh38)'
#  echo -e '-o <-- Absolute path of your output file'
#  echo -e '-v <-- Please select your preferred VEP software version (105 only)'
#  exit 1
#fi

#echo "${@}"
#while getopts ":i:b:o:v:" opt ${@}; do
#        case $opt in
#          i)
#             echo ${OPTARG}
#             infolder=${OPTARG}
#          ;;
#          b)
#             echo ${OPTARG}
#             build=${OPTARG}
#          ;;
#          o)
#             echo ${OPTARG}
#             outfolder=${OPTARG}
#          ;;
#          v)
#             echo ${OPTARG}
#             vep_version=${OPTARG}
#          ;;
#       esac
#done
#
#${vep_bin} --format id -i nome_file.txt --stats_text --force_overwrite  --offline --fasta /shared/resources/hgRef/hg38/Homo_sapiens_assembly38.fasta  --dir_plu
gins ${plugin_path} --o output_file.txt --buffer_size 200000 --fork 12 --plugin CADD,${cadd_whole_genome},${gnomad_indel},${gnomad_v} --custom ${gnomad_vcf},gno
mADv3.1,vcf,exact,0,${gnomad_annot_string} --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --h
gvs --hgvsg --canonical --polyphen b --symbol --vcf --cache --merged --distance 250000 --dir_cache ${maincache}/${vep_version} --assembly ${build} --plugin dbNS
FP,${dbNSFP},${dbNSFP_fields} --plugin Conservation,${gerp_db}
