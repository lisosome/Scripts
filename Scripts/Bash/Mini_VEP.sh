#!/usr/bin/env bash

##Mini VEP script, providing VCFs from a file and adding phyloP and ClinVar information
infile=$1
build=$2
chr=$3
outfolder=$4

## RESOURCES ##
maincache=/shared/resources/VEP_cache
gerp_db=${maincache}/105/compara/gerp_conservation_scores.homo_sapiens.${build}.bw
dbNSFP_fields="Ensembl_transcriptid,Uniprot_acc,VEP_canonical,SIFT4G_score,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE,MetaSVM_score,MetaSVM_rankscore,MetaSVM_pred,MetaLR_score,MetaLR_rankscore,MetaLR_pred,REVEL_score,REVEL_rankscore,ClinPred_score,ClinPred_rankscore,ClinPred_pred,phyloP100way_vertebrate,phyloP100way_vertebrate_rankscore"
gnomad_annot_string="AC,AN,AF,nhomalt,AC-XY,AN-XY,AF-XY,nhomalt-XY,AC-oth,AN-oth,AF-oth,nhomalt-oth,AC-ami,AN-ami,AF-ami,nhomalt-ami,AC-sas,AN-sas,AF-sas,nhomalt-sas,AC-fin,AN-fin,AF-fin,nhomalt-fin,AC-eas,AN-eas,AF-eas,nhomalt-eas,AC-amr,AN-amr,AF-amr,nhomalt-amr,AC-afr,AN-afr,AF-afr,nhomalt-afr,AC-mid,AN-mid,AF-mid,nhomalt-mid,AC-asj,AN-asj,AF-asj,nhomalt-asj,AC-nfe,AN-nfe,AF-nfe,nhomalt-nfe"
#### GRCh37 ####
cadd_v="v1.6"
gnomad_v=${maincache}/105/cadd/${cadd_v}/${build}/gnomad.genomes.r2.1.1.snv.tsv.gz
gnomad_indel=${maincache}/105/cadd/${cadd_v}/${build}/gnomad.genomes.r2.1.1.indel.tsv.gz
cadd_whole_genome=${maincache}/105/cadd/${cadd_v}/${build}/whole_genome_SNVs.tsv.gz
plugin_path=/home/nardone/.conda/envs/vep105/share/ensembl-vep-105.0-0
vep_bin=/home/nardone/.conda/envs/vep105/share/ensembl-vep-105.0-0/vep
gnomad_vcf=/netapp05/analisi_nardone/gnomAD/GRCh37/gnomad.genomes.r2.1.sites.chr${chr}_noVEP.vcf.gz
dbNSFP=/netapp05/analisi_nardone/dbNSFP/${build}/dbNSFP4.3a_grch37_${chr}.gz
ClinVar=/netapp05/analisi_nardone/ClinVar/GRCh37/clinvar.vcf.gz

### Let's set our while loop feeding the input file with the sample code in the first column and the absolute path of VCFs to annotate
source activate vep105
while read line
do
  out_name=$(echo $line | awk '{print $1}')
  VCF=$(echo $line | awk '{print $2}')
  ${vep_bin} -i ${VCF} --stats_text --force_overwrite --format vcf --offline --fasta /shared/resources/gatk4hg19db/Homo_sapiens_assembly19_1000genomes_decoy.fasta --compress_output bgzip --dir_plugins ${plugin_path} --o ${outfolder}/${out_name}.annotated.vcf.gz --buffer_size 200000 --fork 12 --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin CADD,${cadd_whole_genome},${gnomad_indel},${gnomad_v} --plugin Conservation,${gerp_db} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --custom ${ClinVar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --symbol --vcf --cache --merged --dir_cache ${maincache}/105 --assembly ${build}
done < ${infile}
