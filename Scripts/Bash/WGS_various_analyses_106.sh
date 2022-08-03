#!/usr/bin/env bash

#10/11/2017
#script for various analyses on WGS data for the paper, mainly

##USAGE:
#
# WGS_various_analyses.sh <mode> <inputfile> <outputfolder>
#
#MODES:
# SUBSET
# CSQEXTR
# CSQANNOT
# VCFCONCAT
# CSQENRICH
# CADDENRICH
# CADDENRICHAC
# LOFEXTRLAX
# LOFEXTR
# LOFEXTRSTRICT
# SHIFTDET
# SHIFTDETAC
# GET_SHARED
# GENEEXTR
# LOFVCFEXTRLAX

set -e

mode=$1
inputfile=$2
outputfolder=$3
bash_base=/large/___HOME___/burlo/fcrudele/scripts/bash_scripts/

#the mode parameter trigger the reading of other input of different type
case $mode in
	SUBSET )
	#subset the sample in a bootstrap-like way
	#Require:
	# - a sample file with ids
	# - the dimension of the subset
	sample_file=$4
	sample_size=$5
	n_repeat=$6
	storage=$7

	sample_file_name=`basename ${sample_file}`
	output_file_name=`basename ${inputfile}`

	#select a subset of the samples: sort randomly the samples and extract the first sample_size. Repeat it n times.
	m=5G

	for repeat in $(seq 1 ${n_repeat})
	do
		outdir=${outputfolder}/${repeat}
		mkdir -p ${outdir}

		echo "shuf -n ${sample_size} ${sample_file} -o ${outdir}/${sample_file_name%.*}.list" | qsub -N ${mode}_${repeat}_${output_file_name}_samples -o ${outdir}/\$JOB_ID_${mode}_${repeat}_${output_file_name}_samples.log -e ${outdir}/\$JOB_ID_${mode}_${repeat}_${output_file_name}_samples.e -V -l h_vmem=${m} -hold_jid STEP${step_p}_${pop}_${chr}
        echo "bcftools view -S ${outdir}/${sample_file_name%.*}.list ${inputfile} -O z -o ${outdir}/${output_file_name%.vcf*}.vcf.gz" | qsub -N ${mode}_${repeat}_${output_file_name}_vcf -o ${outdir}/\$JOB_ID_${mode}_${repeat}_${output_file_name}_vcf.log -e ${outdir}/\$JOB_ID_${mode}_${repeat}_${output_file_name}_vcf.e -V -l h_vmem=${m} -hold_jid ${mode}_${repeat}_${output_file_name}_samples
		echo "bcftools stats --af-bins 0.01,0.05 -s - -v ${outdir}/${output_file_name%.vcf*}.vcf.gz > ${outdir}/${output_file_name%.vcf*}_stats_subset.tab" | qsub -N ${mode}_${repeat}_${output_file_name}_stats -o ${outdir}/\$JOB_ID_${mode}_${repeat}_${output_file_name}_stats.log -e ${outdir}/\$JOB_ID_${mode}_${repeat}_${output_file_name}_stats.e -V -l h_vmem=${m} -hold_jid ${mode}_${repeat}_${output_file_name}_vcf
		echo "/home/cocca/scripts/bash_scripts/paperONE_pipeline.sh EXTRFREQ ${outdir}/${output_file_name%.vcf*}.vcf.gz" | qsub -N ${mode}_FREQ_${repeat}_${output_file_name}_vcf -o ${outdir}/\$JOB_ID_${mode}_FREQ_${repeat}_${output_file_name}_vcf.log -e ${outdir}/\$JOB_ID_${mode}_FREQ_${repeat}_${output_file_name}_vcf.e -V -l h_vmem=${m} -hold_jid ${mode}_${repeat}_${output_file_name}_vcf

		#clean up a litlle bit, moving stuff into storage
		#default to:
		storage_folder=${storage}/${repeat}
		mkdir -p ${storage_folder}
        echo "/home/cocca/scripts/bash_scripts/paperONE_pipeline.sh CLEANTMP ${outdir}/${output_file_name%.vcf*}.vcf.gz ${outdir}/${output_file_name%.vcf*}.vcf.gz.freq.tab ${storage_folder}/${output_file_name%.vcf*}.vcf.gz.freq.tab"| qsub -N ${mode}_CLEAN_${repeat}_${output_file_name}_vcf -o ${outdir}/\$JOB_ID_${mode}_CLEAN_${repeat}_${output_file_name}_vcf.log -e ${outdir}/\$JOB_ID_${mode}_CLEAN_${repeat}_${output_file_name}_vcf.e -V -l h_vmem=${m} -hold_jid ${mode}_FREQ_${repeat}_${output_file_name}_vcf
		echo "/home/cocca/scripts/bash_scripts/paperONE_pipeline.sh CLEANTMP ${outdir}/${output_file_name%.vcf*}.vcf.gz ${outdir}/${output_file_name%.vcf*}_stats_subset.tab ${storage_folder}/${output_file_name%.vcf*}_stats_subset.tab"| qsub -N ${mode}_CLEAN_${repeat}_${output_file_name}_stats -o ${outdir}/\$JOB_ID_${mode}_CLEAN_${repeat}_${output_file_name}_stats.log -e ${outdir}/\$JOB_ID_${mode}_CLEAN_${repeat}_${output_file_name}_stats.e -V -l h_vmem=${m} -hold_jid ${mode}_FREQ_${repeat}_${output_file_name}_vcf
		# echo "rm ${outdir}/${output_file_name%.vcf*}.vcf.gz;rsync -r -auvzP --exclude=".*" --port=2222 ${outdir}/${output_file_name%.vcf*}.vcf.gz.freq.tab ${storage_folder}/${output_file_name%.vcf*}.vcf.gz.freq.tab;check_and_clean ${outdir}/${output_file_name%.vcf*}.vcf.gz.freq.tab ${storage_folder}/${output_file_name%.vcf*}.vcf.gz.freq.tab GZIP" | qsub -N ${mode}_CLEAN_${repeat}_${output_file_name}_vcf -o ${outdir}/\$JOB_ID_${mode}_CLEAN_${repeat}_${output_file_name}_vcf.log -e ${outdir}/\$JOB_ID_${mode}_CLEAN_${repeat}_${output_file_name}_vcf.e -V -l h_vmem=${m} -hold_jid ${mode}_FREQ_${repeat}_${output_file_name}_vcf

	done

	;;
	CSQEXTR )
	#extract csq categories and frequencies to evaluate enrichment
	chr=$4
	pop=$5
	m=5G

	# 1) extract a list of csq from each chr
	echo "${bash_base}/understand_csq.py ${inputfile} ${outputfolder} ${chr}" | qsub -N ${mode}_${chr}_${pop} -o ${outputfolder}/\$JOB_ID_${mode}_${chr}_${pop}.log -e ${outputfolder}/\$JOB_ID_${mode}_${chr}_${pop}.e -V -l h_vmem=${m}

	#2: merge all lists together
	echo "cat ${outputfolder}/*_consequences.list | sort | uniq > ${outputfolder}/all_consequences.list" |  qsub -N ${mode}_concat_cons_${pop} -o ${outputfolder}/\$JOB_ID_${mode}_concat_cons_${pop}.log -e ${outputfolder}/\$JOB_ID_${mode}_concat_cons_${pop}.e -V -l h_vmem=${m} -hold_jid ${mode}_*_${pop}

	;;

	CSQANNOT)

	infile=${inputfile}
	outfile=$3
	chr=$4
    annot_opt=$5

    cluster=$6
    #add 3 different modes:
    # NEW: brand new annotation of the input file
    # STRIP: stripo old CSQ annotation and reannotate
    # REGION: select a region for annotation
    #define vep cache basefolder
    # maincache=/tmp/VEP_cache
    # main_sw=z/home/cocca/softwares
    # ref_fasta=${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
    vep_version=$7
    assembly_build=$8
    ref_fasta=$9
    case ${cluster} in
        CBM )
            vep_bin=/home/cocca/softwares/VEP/ensembl-vep/vep
            maincache=/home/cocca/resources/VEP_cache
        ;;
        BURLO)
            maincache=/shared/resources/VEP_cache
            case ${vep_version} in
                90)
                    vep_bin=/shared/software/VEP/ensembl-vep/vep
                    loftee_path=/shared/softwares/VEP/loftee
                ;;
                94)
                    vep_bin=/share/apps/bio/miniconda2/bin/vep
                    loftee_path=/share/apps/bio/miniconda2/share/ensembl-vep-94.5-0
                ;;
                100)
                    #we need to activate the vep env
                    source activate vep
                    vep_bin=/shared/software/VEP/ensembl-vep-100/ensembl-vep/vep
                    loftee_path=/shared/software/VEP/loftee
                    # dbNSFP=/netapp02/data/resources/dbNSFP/dbNSFP4.1a.txt.gz
                    dbNSFP=/netapp02/data/resources/dbNSFP/v4.1a/dbNSFP4.1a_variant.chr${chr}.gz
                    dbNSFP_fields="Ensembl_transcriptid,Uniprot_acc,VEP_canonical,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE"
                    # gnomad_vcf=/fast/burlo/fcrudele/resources/gnomAD/v3.1/gnomad.genomes.v3.1.sites.chr${chr}.vcf.bgz
                ;;
            esac
        ;;
        ORFEO)
            maincache=/fast/burlo/fcrudele/resources/VEP_cache
            gnomad_vcf=/fast/burlo/fcrudele/resources/gnomAD/v3.1/gnomad.genomes.v3.1.sites.chr${chr}.vcf.bgz

            case ${vep_version} in
                90)
                    vep_bin=/shared/software/VEP/ensembl-vep/vep
                    loftee_path=/large/___HOME___/burlo/cocca/software/loftee
                ;;
                94)
                    vep_bin=/share/apps/bio/miniconda2/bin/vep
                    loftee_path=/large/___HOME___/burlo/cocca/software/loftee
                ;;
                106)
                    #we need to activate the vep env
                    eval "$(conda shell.bash hook)"
                    conda activate VEPP
                    vep_bin=/large/___HOME___/burlo/nardone/.conda/envs/VEPP/bin/vep
                    loftee_path=/large/___HOME___/burlo/nardone/.conda/envs/VEPP/share/ensembl-vep-106.1-0
                    dbNSFP=/fast/burlo/fcrudele/resources/dbNSFP/v4.1a/dbNSFP4.1a_variant.chr${chr}.gz
                    dbNSFP_fields="Ensembl_transcriptid,Uniprot_acc,VEP_canonical,Polyphen2_HDIV_score,Polyphen2_HDIV_rankscore,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_rankscore,Polyphen2_HVAR_pred,LRT_score,LRT_converted_rankscore,LRT_pred,LRT_Omega,MutationTaster_score,MutationTaster_converted_rankscore,MutationTaster_pred,MutationTaster_model,MutationTaster_AAE"
                ;;
            esac
        ;;
    esac
    gerp_db=${maincache}/${vep_version}/compara/gerp_conservation_scores.homo_sapiens.${assembly_build}.bw


    echo -e "infile: ${infile}
    outfile: ${outfile}
    chr: ${chr}
    annot_opt: ${annot_opt}
    cluster: ${cluster}
    vep_version: ${vep_version}
    assembly_build: ${assembly_build}
    ref_fasta: ${ref_fasta}
    vep binary: ${vep_bin}
    gerp_db: ${gerp_db}"

maxsize=1000000000000000000000000000000000000000000000000000000
    case $annot_opt in
        NEW)
            # ${main_sw}/ensembl-vep/vep --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/90/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/90 --assembly GRCh37
            # ${main_sw}/ensembl-vep/vep --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/90/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --plugin CADD,${maincache}/90/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/90 --assembly GRCh37
            # ${vep_bin} --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
            # ${vep_bin} --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/v1.4/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/InDels.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/gnomad.genomes.r2.0.1.sites.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
            case ${assembly_build} in
                GRCh38 )
                    gnomad_v=gnomad.genomes.r3.0.snv.tsv.gz
                    indels=gnomad.genomes.r3.0.indel.tsv.gz
                    cadd_v="v1.6"
                    gnomad_annot_string="AC,AN,AF,nhomalt,AC-XY,AN-XY,AF-XY,nhomalt-XY,AC-oth,AN-oth,AF-oth,nhomalt-oth,AC-ami,AN-ami,AF-ami,nhomalt-ami,AC-sas,AN-sas,AF-sas,nhomalt-sas,AC-fin,AN-fin,AF-fin,nhomalt-fin,AC-eas,AN-eas,AF-eas,nhomalt-eas,AC-amr,AN-amr,AF-amr,nhomalt-amr,AC-afr,AN-afr,AF-afr,nhomalt-afr,AC-mid,AN-mid,AF-mid,nhomalt-mid,AC-asj,AN-asj,AF-asj,nhomalt-asj,AC-nfe,AN-nfe,AF-nfe,nhomalt-nfe"

                    #removed gnomad custom annotation to speed up the annotation on all chromosomes
                    if [ -z ${gnomad_vcf+x} ]
                    then
                        echo "gnomad vcf annotations are unset"
                        ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --buffer_size 200000 --fork 12 --fasta ${ref_fasta} --max_sv_size ${maxsize} --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
                    else
                        echo "gnomad vcf annotations are set"
                        ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --buffer_size 200000 --fork 12 --fasta ${ref_fasta} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af_gnomad --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
                    fi
                    # ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --buffer_size 200000 --fork 24 --fasta ${ref_fasta} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin dbNSFP,${dbNSFP},${dbNSFP_fields} --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
                    # ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --buffer_size 200000 --fork 24 --fasta ${ref_fasta} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
                    # ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
                    ;;
                GRCh37 )
                    gnomad_v=gnomad.genomes.r2.0.1.sites.tsv.gz
                    cadd_v="v1.4"
                    indels=InDels.tsv.gz


                    # ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
                    ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --buffer_size 200000 --fork 12 --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
                    ;;
            esac
            echo -e "gnomad_version: ${gnomad_v}
                     cadd_version: ${cadd_v}"

            echo -e "Command issued:\n
            ${vep_bin} -i ${infile} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --fasta ${ref_fasta} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}"


            # ${main_sw}/ensembl-vep/vep --i ${infile} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf.sql,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
        ;;
        STRIP)
            # bcftools annotate -x INFO/CSQ ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
            # bcftools annotate -x INFO/CSQ ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
            # bcftools annotate -x INFO/CSQ ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
            # bcftools annotate -x INFO/CSQ ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/v1.4/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/InDels.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/gnomad.genomes.r2.0.1.sites.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
            # bcftools annotate -x INFO/CSQ ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/v1.4/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/InDels.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/gnomad.genomes.r2.0.1.sites.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
            case ${assembly_build} in
                GRCh38 )
                    gnomad_v=gnomad.genomes.r3.0.snv.tsv.gz
                    indels=gnomad.genomes.r3.0.indel.tsv.gz
                    cadd_v="v1.6"
                    gnomad_annot_string="AC,AN,AF,nhomalt,AC-XY,AN-XY,AF-XY,nhomalt-XY,AC-oth,AN-oth,AF-oth,nhomalt-oth,AC-ami,AN-ami,AF-ami,nhomalt-ami,AC-sas,AN-sas,AF-sas,nhomalt-sas,AC-fin,AN-fin,AF-fin,nhomalt-fin,AC-eas,AN-eas,AF-eas,nhomalt-eas,AC-amr,AN-amr,AF-amr,nhomalt-amr,AC-afr,AN-afr,AF-afr,nhomalt-afr,AC-mid,AN-mid,AF-mid,nhomalt-mid,AC-asj,AN-asj,AF-asj,nhomalt-asj,AC-nfe,AN-nfe,AF-nfe,nhomalt-nfe"
                    ;;
                GRCh37 )
                    gnomad_v=gnomad.genomes.r2.0.1.sites.tsv.gz
                    cadd_v="v1.4"
                    indels=InDels.tsv.gz
                    ;;
            esac
            echo -e "gnomad_version: ${gnomad_v}
                     cadd_version: ${cadd_v}"

            echo -e "Command issued:\n
            bcftools annotate -x INFO/CSQ ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --fasta ${ref_fasta} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}"

            # bcftools annotate -x INFO/CSQ ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
            bcftools annotate -x INFO/CSQ ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --buffer_size 200000 --fork 16 --fasta ${ref_fasta} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
            # bcftools annotate -x INFO/CSQ ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/${vep_version}/homo_sapiens/${vep_version}_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf.sql,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
        ;;
        REGION)
            case ${assembly_build} in
                GRCh38 )
                    gnomad_v=gnomad.genomes.r3.0.snv.tsv.gz
                    indels=gnomad.genomes.r3.0.indel.tsv.gz
                    cadd_v="v1.6"
                    gnomad_annot_string="AC,AN,AF,nhomalt,AC-XY,AN-XY,AF-XY,nhomalt-XY,AC-oth,AN-oth,AF-oth,nhomalt-oth,AC-ami,AN-ami,AF-ami,nhomalt-ami,AC-sas,AN-sas,AF-sas,nhomalt-sas,AC-fin,AN-fin,AF-fin,nhomalt-fin,AC-eas,AN-eas,AF-eas,nhomalt-eas,AC-amr,AN-amr,AF-amr,nhomalt-amr,AC-afr,AN-afr,AF-afr,nhomalt-afr,AC-mid,AN-mid,AF-mid,nhomalt-mid,AC-asj,AN-asj,AF-asj,nhomalt-asj,AC-nfe,AN-nfe,AF-nfe,nhomalt-nfe"
                    ;;
                GRCh37 )
                    gnomad_v=gnomad.genomes.r2.0.1.sites.tsv.gz
                    cadd_v="v1.4"
                    indels=InDels.tsv.gz
                    ;;
            esac
            echo -e "gnomad_version: ${gnomad_v}
                     cadd_version: ${cadd_v}"

            start=${10}
            end=${11}
            # bcftools view -r ${chr}:${start}-${end} ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
            # bcftools view -r ${chr}:${start}-${end} ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,loftee_path:${main_sw}/loftee --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly GRCh37
            # bcftools view -r ${chr}:${start}-${end} ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd13/1000G_phase3.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
            # bcftools view -r ${chr}:${start}-${end} ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/v1.4/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/InDels.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/gnomad.genomes.r2.0.1.sites.tsv.gz --symbol --vcf --cache --merged --dir ${maincache}/${vep_version} --assembly ${assembly_build}
            #need to get the correct naming for the current chromosome:
            chr_r=$(bcftools view -H ${infile}|head -1|cut -f 1)
            # bcftools view -r ${chr}:${start}-${end} ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/v1.4/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/InDels.tsv.gz,${maincache}/${vep_version}/cadd/v1.4/gnomad.genomes.r2.0.1.sites.tsv.gz --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
            bcftools view -r ${chr_r}:${start}-${end} ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --buffer_size 200000 --fork 16 --fasta ${ref_fasta} --custom ${gnomad_vcf},gnomADv3.1,vcf,exact,0,${gnomad_annot_string} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --xref_refseq --hgvs --hgvsg --canonical --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
            # bcftools view -r ${chr}:${start}-${end} ${infile} | ${vep_bin} --stats_text --force_overwrite --format vcf --offline --compress_output bgzip --o ${outfile} --fasta ${ref_fasta} --af --af_1kg --af_esp --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/${vep_version}/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/${vep_version}/loftee_data/phylocsf_gerp.sql,loftee_path:${loftee_path} --plugin CADD,${maincache}/${vep_version}/cadd/${cadd_v}/whole_genome_SNVs.tsv.gz,${maincache}/${vep_version}/cadd/${cadd_v}/${indels},${maincache}/${vep_version}/cadd/${cadd_v}/${gnomad_v} --plugin Conservation,${gerp_db} --symbol --vcf --cache --merged --dir_cache ${maincache}/${vep_version} --assembly ${assembly_build}
            # bcftools view -r ${chr}:${start}-${end} ${infile} | ${main_sw}/ensembl-vep/vep --stats_text --force_overwrite --format vcf --offline --o ${outfile} --fasta ${maincache}/90/homo_sapiens/90_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz --af --af_1kg --af_esp --af_gnomad --variant_class --regulatory --ccds --protein --uniprot --sift b --polyphen b --plugin LoF,human_ancestor_fa:${maincache}/90/loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15,conservation_file:${maincache}/90/loftee_data/phylocsf.sql,loftee_path:${main_sw}/loftee --symbol --vcf --cache --merged --dir ${maincache}/90 --assembly GRCh37
        ;;
    esac

    # annotate vcf files with consequences
    echo "index output file"
	# bgzip ${outfile}
	tabix -p vcf ${outfile}
	;;
	VCFCONCAT)
	#concat vcf files
	file_list=${inputfile}
	outfile=$3

	bcftools concat -f ${file_list} -O z -o ${outfile}.unsorted.vcf.gz
	bcftools sort ${outfile}.unsorted.vcf.gz -O z -o ${outfile}
	tabix -p vcf ${outfile}
	;;
	CSQENRICH)

	chr=$4
	pop=$5
    file_type=$6
	# inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
	# outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
	input_name=`basename ${inputfile}`

	mkdir -p ${outputfolder}
	#extract info from VCF files to count variants to use a maf splitter style function
	# (echo "CHROM POS ID REF ALT AC AN IMP2 VQSLOD AF MAF MINOR";
    #add a little bit of smartnessness: if the extension of the file is vcf.gz or vcf, it means we have to generate the table, otherwise, we assume the table is already the inputfile!
    # We want to manage vcf, vcf.gz

    case $file_type in
        VCF )
            (bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\n' ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        VCFGZ )
            (bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\n' ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        TAB )
            # it means the tab file is the input file...
            (zcat ${inputfile} | cut -f -8 -d " " | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
            # (zcat ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
    esac


    #generate the baseline freq table
    # maf_bins="0,5,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50"
    maf_bins="0,4,0.02,0.05,0.10,0.20,0.30,0.40,0.50"
    # maf_bins="0,4,0.02,0.05,0.10,0.50"

    (echo "CHROM POS ID REF ALT AC AN CSQ AF MAF MIN COHORT MAF BIN";(zgrep -v CHROM ${outputfolder}/${input_name}.maf_csq.tab.gz | awk -v cohort=$pop '{if($9 <= 0.5 ) print $0,toupper(cohort),$9; else print $0,toupper(cohort),1-$9}'|awk -v bins=$maf_bins '
    {n=split(bins,mafs,",");}{
    for (i=1;i<=n;i++){
        if (mafs[i]==0){
            if($(NF) <= mafs[i]){
            print $0,mafs[i]
            }
        }else{
            if (mafs[i]==4){
                if( ($7-$6) <= mafs[i] && ($7-$6) > mafs[i-1]){
                    print $0,mafs[i]
                }
            }else{
                if(mafs[i-1]==4){
                    if($(NF) <= mafs[i] && ($7-$6) > mafs[i-1]){
                        print $0,mafs[i]
                    }
                }else{
                    if($(NF) <= mafs[i] && $(NF) > mafs[i-1]){
                        print $0,mafs[i]
                    }
                }
            }
        }
    }
    }'))| gzip -c > ${outputfolder}/${input_name}.maf_csq_bin.tab.gz



    #extract counts for each freq bin
    zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk '{print $(NF)}'| sort -g| uniq -c | fgrep -v BIN | awk '{print $1,$2}' > ${outputfolder}/${input_name}.maf_csq_all_bin_resume.tab

    #define lof variants
    declare -a lof_group
    lof_group=(frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_gained stop_lost)
    mkdir -p ${outputfolder}/lof_cons
    # for lof in ${lof_group[@]}
    for (( il = 0; il < ${#lof_group[@]}; il++ ))
    do
        lof=${lof_group[il]}
        case $lof in
            frameshift_variant )
                zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            splice_acceptor_variant )
                zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            splice_donor_variant )
                zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            start_lost )
                zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            stop_gained )
                zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            stop_lost )
                zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" ' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
        esac

        # zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | fgrep $lof | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    done

    for maf_bin in 0 4 0.02 0.05 0.10 0.20 0.30 0.40 0.50
    do

    zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | awk -v mafbin=${maf_bin} '{if($(NF)==mafbin) print $(NF)}'| sort -g| uniq -c | awk '{print $1,$2}'> ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab

    done

    #stratify by consequence category
    declare -a consequences
    consequences=(synonymous_variant missense_variant)
    # for cond in ${consequences[@]}
    for (( i = 0; i < ${#consequences[@]}; i++ ))
    do
        # for this set, we need to remove all the data from those set with higher pathogenicity:
        # so, if we want to count the synonymous variants we need to check if a variant is counted as a LOF, and we need to remove it from out set,
        #  than we need to remove all variants that are also annotated somewhere as missense_variant the same we have to do with the other category

        cond=${consequences[i]}
        if [[ $i -eq 0 ]]; then
            #statement
            n_cond=${consequences[i+1]}
        else
            n_cond=${consequences[i-1]}
        fi
        zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v conseq=${cond} -v n_conseq=${n_cond} '$8~conseq && $8!~n_conseq && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/${input_name}.maf_${cond}_bin.tab.gz
        zcat ${outputfolder}/${input_name}.maf_${cond}_bin.tab.gz | awk '{print $(NF)}'| sort -g| uniq -c | awk '{print $1,$2}'> ${outputfolder}/${input_name}.maf_${cond}_bin_resume.tab
    done

    ;;
    CADDENRICH)

    chr=$4
    pop=$5
    file_type=$6
    # inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
    # outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
    input_name=`basename ${inputfile}`

    mkdir -p ${outputfolder}
    #extract info from VCF files to count variants to use a maf splitter style function
    # (echo "CHROM POS ID REF ALT AC AN IMP2 VQSLOD AF MAF MINOR";
    #add a little bit of smartnessness: if the extension of the file is vcf.gz or vcf, it means we have to generate the table, otherwise, we assume the table is already the inputfile!
    # We want to manage vcf, vcf.gz

    case $file_type in
        VCF )
            (bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\t%INFO/CADD_PHRED\n' ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        VCFGZ )
            (bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\t%INFO/CADD_PHRED\n' ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        TAB )
            # it means the tab file is the input file...
            (zcat ${inputfile} | cut -f -9 -d " " | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
            # (zcat ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        VEP90 )
            # it means the tab file is the input file...
            ~/scripts/bash_scripts/vep_parser.py --input ${inputfile} --cadd --out ${outputfolder}/${input_name}.tab.gz
            #here we need to add a NA in the CSQ column, for consistency with the other file types
            (zcat ${outputfolder}/${input_name}.tab.gz | awk '{print $1,$2,$3,$4,$5,$6,$7,"NA",$8}' | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
            # (zcat ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
    esac


    #generate the baseline freq table
    # maf_bins="0,5,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50"
    maf_bins="0,4,0.02,0.05,0.10,0.20,0.30,0.40,0.50"

    # maf_bins="0,4,0.02,0.05,0.10,0.50"
    cadd_bins="0,5,15,20"

    (echo "CHROM POS ID REF ALT AC AN CSQ CADD_PHRED AF MAF MIN COHORT MAF MAF_BIN CADD_BIN";(zgrep -v CHROM ${outputfolder}/${input_name}.maf_csq.tab.gz | awk -v cohort=$pop '{if($10 <= 0.5 ) print $0,toupper(cohort),$10; else print $0,toupper(cohort),1-$10}'|awk -v bins=${maf_bins} -v cadds=${cadd_bins} '
    {n=split(bins,mafs,",");c=split(cadds,cadd,",")}{
    for (i=1;i<=n;i++){
        for (j=1;j<=c;j++){
            if (cadd[j]==20){
                if ( $9 >= cadd[j] ){
                    if (mafs[i]==0){
                        if($(NF) <= mafs[i]){
                            print $0,mafs[i],cadd[j]
                        }
                    }else{
                        if (mafs[i]==4){
                            if(($7-$6) <= mafs[i] && ($7-$6) > mafs[i-1]){
                                print $0,mafs[i],cadd[j]
                            } else if ($6 <= mafs[i] && $6 > mafs[i-1]) {
                                print $0,mafs[i],cadd[j]
                            }
                        }else{
                            if(mafs[i-1]==4){
                                if($(NF) <= mafs[i] && ($7-$6) > mafs[i-1]){
                                    print $0,mafs[i],cadd[j]
                                } else if ($NF <= mafs[i] && $6 > mafs[i-1]) {
                                    print $0,mafs[i],cadd[j]
                                }
                            }else{
                                if($(NF) <= mafs[i] && $(NF) > mafs[i-1]){
                                    print $0,mafs[i],cadd[j]
                                }
                            }
                        }
                    }
                }
            }else{
                if($9 >= cadd[j] && $9 < cadd[j+1]) {
                    if (mafs[i]==0){
                        if($(NF) <= mafs[i]){
                        print $0,mafs[i],cadd[j]"_"cadd[j+1]
                        }
                    }else{
                        if (mafs[i]==4){
                            if(($7-$6) <= mafs[i] && ($7-$6) > mafs[i-1]){
                                print $0,mafs[i],cadd[j]"_"cadd[j+1]
                            } else if ($6 <= mafs[i] && $6 > mafs[i-1]) {
                                print $0,mafs[i],cadd[j]"_"cadd[j+1]
                            }
                        }else{
                            if(mafs[i-1]==4){
                                if($(NF) <= mafs[i] && ($7-$6) > mafs[i-1]){
                                    print $0,mafs[i],cadd[j]"_"cadd[j+1]
                                } else if ($NF <= mafs[i] && $6 > mafs[i-1]) {
                                    print $0,mafs[i],cadd[j]"_"cadd[j+1]
                                }
                            }else{
                                if($(NF) <= mafs[i] && $(NF) > mafs[i-1]){
                                    print $0,mafs[i],cadd[j]"_"cadd[j+1]
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    }'))| gzip -c > ${outputfolder}/${input_name}.maf_csq_bin.tab.gz



    #extract counts for each freq bin
    zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk '{print $(NF-1)}'| sort -g| uniq -c | fgrep -v BIN | awk '{print $1,$2}' > ${outputfolder}/${input_name}.maf_csq_all_bin_resume.tab


    #extract counts for each freq bin and for each CADD_PHRED range

    #define lof variants
    declare -a cadd_group
    cadd_group=(0_5 5_15 15_20 20)
    mkdir -p ${outputfolder}/cadd_cons
    # for lof in ${lof_group[@]}
    for (( il = 0; il < ${#cadd_group[@]}; il++ ))
    do
        cadd=${cadd_group[il]}
        zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v cadd_range=${cadd} '$(NF)==cadd_range' | gzip -c >  ${outputfolder}/cadd_cons/${input_name}.cadd_${cadd}_bin.tab.gz
    done

    #resume counts for each freq bin and for each CADD_PHRED range
    zcat ${outputfolder}/cadd_cons/${input_name}.cadd_*_bin.tab.gz | awk '{print $(NF-1),$(NF)}'| sort -g| uniq -c | awk '{print $1,$2,$3}'> ${outputfolder}/cadd_cons/${input_name}.cadd_all_maf_bin_resume.tab

    ;;
    CADDENRICHAC)

    chr=$4
    pop=$5
    file_type=$6
    # inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
    # outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
    input_name=`basename ${inputfile}`

    mkdir -p ${outputfolder}
    #extract info from VCF files to count variants to use a maf splitter style function
    # (echo "CHROM POS ID REF ALT AC AN IMP2 VQSLOD AF MAF MINOR";
    #add a little bit of smartnessness: if the extension of the file is vcf.gz or vcf, it means we have to generate the table, otherwise, we assume the table is already the inputfile!
    # We want to manage vcf, vcf.gz

    case $file_type in
        VCF )
            (bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\t%INFO/CADD_PHRED\n' ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        VCFGZ )
            (bcftools query -f'%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/CSQ\t%INFO/CADD_PHRED\n' ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        TAB )
            # it means the tab file is the input file...
            (zcat ${inputfile} | cut -f -9 -d " " | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
            # (zcat ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
        VEP90 )
            # it means the tab file is the input file...
            ~/scripts/bash_scripts/vep_parser.py --input ${inputfile} --cadd --out ${outputfolder}/${input_name}.tab.gz
            #here we need to add a NA in the CSQ column, for consistency with the other file types
            (zcat ${outputfolder}/${input_name}.tab.gz | awk '{print $1,$2,$3,$4,$5,$6,$7,"NA",$8}' | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
            # (zcat ${inputfile} | fgrep -v CHROM | awk '{if($5 !~ ",") print $0}' | awk '{if($7 != 0)print $0,$6/$7;else print $0,"NA"}' | awk '{if($NF != "NA") {if ($NF <= 0.5) print $0,$NF,$5;else print $0,(1-$NF),$4}else{print $0,"NA","NA"}}')| tr "\t" " "| uniq | gzip -c > ${outputfolder}/${input_name}.maf_csq.tab.gz
        ;;
    esac


    #generate the baseline freq table
    # maf_bins="0,5,0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50"
    maf_bins="0,2,5,10,10"

    # maf_bins="0,4,0.02,0.05,0.10,0.50"
    cadd_bins="0,5,15,20"

    (echo "CHROM POS ID REF ALT AC AN CSQ CADD_PHRED AF MAF MIN COHORT MAF MAF_BIN CADD_BIN";(zgrep -v CHROM ${outputfolder}/${input_name}.maf_csq.tab.gz | awk -v cohort=$pop '{if($10 <= 0.5 ) print $0,toupper(cohort),$10; else print $0,toupper(cohort),1-$10}'|awk -v bins=${maf_bins} -v cadds=${cadd_bins} '
    {n=split(bins,mafs,",");c=split(cadds,cadd,",")}{
    for (i=1;i<=n;i++){
        for (j=1;j<=c;j++){
            if (cadd[j]==20){
                if ( $9 >= cadd[j] ){
                    if (i==1){
                        if($6 <= mafs[i]){
                            print $0,mafs[i],cadd[j]
                        }
                    }else if (i==n) {
                        if($6 > mafs[i]) {
                            print $0,mafs[i]"_L",cadd[j]
                        }
                    }else {
                        if($6 <= mafs[i] && $6 > mafs[i-1]) {
                            print $0,mafs[i],cadd[j]
                        }
                    }
                }
            }else{
                if($9 >= cadd[j] && $9 < cadd[j+1]) {
                    if (i==1){
                        if($6 <= mafs[i]){
                            print $0,mafs[i],cadd[j]"_"cadd[j+1]
                        }
                    }else if (i==n) {
                        if($6 > mafs[i]) {
                            print $0,mafs[i]"_L",cadd[j]"_"cadd[j+1]
                        }
                    }else {
                        if($6 <= mafs[i] && $6 > mafs[i-1]) {
                            print $0,mafs[i],cadd[j]"_"cadd[j+1]
                        }
                    }
                }
            }
        }
    }
    }'))| gzip -c > ${outputfolder}/${input_name}.maf_csq_bin.tab.gz



    #extract counts for each freq bin
    zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk '{print $(NF-1)}'| sort -g| uniq -c | fgrep -v BIN | awk '{print $1,$2}' > ${outputfolder}/${input_name}.maf_csq_all_bin_resume.tab


    #extract counts for each freq bin and for each CADD_PHRED range

    #define lof variants
    declare -a cadd_group
    cadd_group=(0_5 5_15 15_20 20)
    mkdir -p ${outputfolder}/cadd_cons
    # for lof in ${lof_group[@]}
    for (( il = 0; il < ${#cadd_group[@]}; il++ ))
    do
        cadd=${cadd_group[il]}
        zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | awk -v cadd_range=${cadd} '$(NF)==cadd_range' | gzip -c >  ${outputfolder}/cadd_cons/${input_name}.cadd_${cadd}_bin.tab.gz
    done

    #resume counts for each freq bin and for each CADD_PHRED range
    zcat ${outputfolder}/cadd_cons/${input_name}.cadd_*_bin.tab.gz | awk '{print $(NF-1),$(NF)}'| sort -g| uniq -c | awk '{print $1,$2,$3}'> ${outputfolder}/cadd_cons/${input_name}.cadd_all_maf_bin_resume.tab
    ;;
    LOFEXTRLAX )
    # mode=$1
    # inputfile=$2
    # outputfolder=$3

    chr=$4
    pop=$5
    ref_fasta=$6
    # file_type=$6
    # inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
    # outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
    input_name=`basename ${inputfile}`

    #define lof variants
    declare -a lof_group
    lof_group=(frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_gained stop_lost)
    mkdir -p ${outputfolder}/lof_cons

    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -f ${ref_fasta} -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%INFO/CADD_PHRED\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz

    # for lof in ${lof_group[@]}
    for (( il = 0; il < ${#lof_group[@]}; il++ ))
    do
        lof=${lof_group[il]}
        # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
        # zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
        zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
        # zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | fgrep $lof | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    done
    #we need also to create a merged version of the data extracted
    zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | sort -k1,1 -k2,2 | uniq | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab.gz

    #clean splitted file
    rm ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    ;;
    LOFEXTR )
    # mode=$1
    # inputfile=$2
    # outputfolder=$3

    chr=$4
    pop=$5
    # file_type=$6
    # inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
    # outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
    mkdir -p ${outputfolder}/lof_cons
    input_name=`basename ${inputfile}`

    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -f ${ref_fasta} -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%INFO/CADD_PHRED\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz

    #define lof variants
    declare -a lof_group
    lof_group=(frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_gained stop_lost)
    # for lof in ${lof_group[@]}
    for (( il = 0; il < ${#lof_group[@]}; il++ ))
    do
        lof=${lof_group[il]}
        case $lof in
            frameshift_variant )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            splice_acceptor_variant )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            splice_donor_variant )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            start_lost )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            stop_gained )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            stop_lost )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" ' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" ' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
        esac

        # zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | fgrep $lof | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    done

    zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | awk '{print $(NF)}'| sort -g| uniq -c | awk '{print $1,$2}'> ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab
    #we need also to create a merged version of the data extracted
    zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | sort -k1,1 -k2,2 | uniq | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab.gz

    #stratify by consequence category
    declare -a consequences
    consequences=(synonymous_variant missense_variant)
    # for cond in ${consequences[@]}
    for (( i = 0; i < ${#consequences[@]}; i++ ))
    do
        # for this set, we need to remove all the data from those set with higher pathogenicity:
        # so, if we want to count the synonymous variants we need to check if a variant is counted as a LOF, and we need to remove it from out set,
        #  than we need to remove all variants that are also annotated somewhere as missense_variant the same we have to do with the other category

        cond=${consequences[i]}
        if [[ $i -eq 0 ]]; then
            #statement
            n_cond=${consequences[i+1]}
        else
            n_cond=${consequences[i-1]}
        fi
        zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${cond} -v n_conseq=${n_cond} '$8~conseq && $8!~n_conseq && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/${input_name}.maf_${cond}_bin.tab.gz
        zcat ${outputfolder}/${input_name}.maf_${cond}_bin.tab.gz | awk '{print $(NF)}'| sort -g| uniq -c | awk '{print $1,$2}'> ${outputfolder}/${input_name}.maf_${cond}_bin_resume.tab
    done

    #clean splitted file
    rm ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz

    ;;
    LOFEXTRSTRICT )
    # mode=$1
    # inputfile=$2
    # outputfolder=$3

    chr=$4
    pop=$5
    # file_type=$6
    # inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
    # outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
    mkdir -p ${outputfolder}/lof_cons
    input_name=`basename ${inputfile}`

    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -f ${ref_fasta} -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%INFO/CADD_PHRED\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz

    #define lof variants
    declare -a lof_group
    lof_group=(frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_gained stop_lost)
    # for lof in ${lof_group[@]}
    for (( il = 0; il < ${#lof_group[@]}; il++ ))
    do
        lof=${lof_group[il]}
        case $lof in
            frameshift_variant )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"intron_variant" && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            splice_acceptor_variant )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"intron_variant" && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            splice_donor_variant )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"intron_variant" && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            start_lost )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"intron_variant" && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            stop_gained )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"intron_variant" && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
            stop_lost )
                # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" ' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
                zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"intron_variant" && $8!~"synonymous_variant" && $8!~"missense_variant" && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" ' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
            ;;
        esac

        # zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | fgrep $lof | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    done

    zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | awk '{print $(NF)}'| sort -g| uniq -c | awk '{print $1,$2}'> ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab
    #we need also to create a merged version of the data extracted
    zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | sort -k1,1 -k2,2 | uniq | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab.gz

    #stratify by consequence category
    declare -a consequences
    consequences=(synonymous_variant missense_variant)
    # for cond in ${consequences[@]}
    for (( i = 0; i < ${#consequences[@]}; i++ ))
    do
        # for this set, we need to remove all the data from those set with higher pathogenicity:
        # so, if we want to count the synonymous variants we need to check if a variant is counted as a LOF, and we need to remove it from out set,
        #  than we need to remove all variants that are also annotated somewhere as missense_variant the same we have to do with the other category

        cond=${consequences[i]}
        if [[ $i -eq 0 ]]; then
            #statement
            n_cond=${consequences[i+1]}
        else
            n_cond=${consequences[i-1]}
        fi
        zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${cond} -v n_conseq=${n_cond} '$8~conseq && $8!~n_conseq && $8!~"frameshift_variant" && $8!~"splice_acceptor_variant" && $8!~"splice_donor_variant" && $8!~"start_lost" && $8!~"stop_gained" && $8!~"stop_lost"' | gzip -c >  ${outputfolder}/${input_name}.maf_${cond}_bin.tab.gz
        zcat ${outputfolder}/${input_name}.maf_${cond}_bin.tab.gz | awk '{print $(NF)}'| sort -g| uniq -c | awk '{print $1,$2}'> ${outputfolder}/${input_name}.maf_${cond}_bin_resume.tab
    done

    #clean splitted file
    rm ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz

    # mode=$1
    # inputfile=$2
    # outputfolder=$3

    # chr=$4
    # pop=$5
    # ref_fasta=$6
    # # file_type=$6
    # # inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
    # # outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
    # input_name=`basename ${inputfile}`

    # #define lof variants
    # declare -a lof_group
    # lof_group=(frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_gained stop_lost)
    # mkdir -p ${outputfolder}/lof_cons

    # # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    # # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -f ${ref_fasta} -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    # # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%INFO/CADD_PHRED\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz

    # # for lof in ${lof_group[@]}
    # for (( il = 0; il < ${#lof_group[@]}; il++ ))
    # do
    #     lof=${lof_group[il]}
    #     # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    #     # zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    #     zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    #     # zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | fgrep $lof | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    # done
    # #we need also to create a merged version of the data extracted
    # zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | sort -k1,1 -k2,2 | uniq | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab.gz

    # #clean splitted file
    # rm ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz


    ;;
    SHIFTDET)
    # mode=$1
    # inputfile=$2
    # outputfolder=$3
    output_file=$3
    other_file=$4
    maf=$5

    zcat ${inputfile} | awk 'FNR==NR {a[$1"_"$2"_"$4"_"$5]=$0;next} $1"_"$2"_"$4"_"$5 in a {print a[$1"_"$2"_"$4"_"$5], $(NF-3),$(NF-2), $(NF-1),$(NF)}' <( zcat ${other_file} | awk -v maf_thr=${maf} '$(NF-2)<=maf_thr' ) - | awk '{ if($(NF-1) != $(NF-5)) print $0,$(NF-6)-$(NF-2),"SHIFTED";else print $0,$(NF-6)-$(NF-2),"SAME"}'| gzip -c > ${output_file}

    ;;
    SHIFTDETAC)
    # mode=$1
    # inputfile=$2
    # outputfolder=$3
    output_file=$3
    other_file=$4
    # maf=$5

    zcat ${inputfile} | awk 'FNR==NR {a[$1"_"$2"_"$4"_"$5]=$0;next} $1"_"$2"_"$4"_"$5 in a {print a[$1"_"$2"_"$4"_"$5], $(NF-3),$(NF-2), $(NF-1),$(NF)}' <( zcat ${other_file} ) - | awk '{ if($(NF-1) != $(NF-5)) print $0,$(NF-6)-$(NF-2),"SHIFTED";else print $0,$(NF-6)-$(NF-2),"SAME"}'| gzip -c > ${output_file}

    ;;
    GET_SHARED )
        # mode=$1
        # inputfile=$2
        # outputfolder=$3
        inputfile2=$4
        pop1=$5
        pop2=$6
        chr=$7
        format=$8
        pop1_folder=`dirname ${inputfile}`
        pop2_folder=`dirname ${inputfile2}`

        common_out_path_pop1=${outputfolder}/${pop1}
        common_out_path_pop2=${outputfolder}/${pop2}

        mkdir -p ${common_out_path_pop1}
        mkdir -p ${common_out_path_pop2}
        mkdir -p ${pop1_folder}/${pop2}
        mkdir -p ${pop2_folder}/${pop1}

        echo "Working on ${chr} in ${pop1} and ${pop2}..."

        #extract relevant columns: Input files have to be in a table format already or a vcf file
        case ${format} in
            VCF )
                bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" ${inputfile} | tr "\t" "_" > ${pop1_folder}/${pop2}/${chr}.${pop1}.tab
                bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" ${inputfile2} | tr "\t" "_" > ${pop2_folder}/${pop1}/${chr}.${pop2}.tab
            ;;
            TAB )
                zcat ${inputfile} | tail -n+2 | awk '{print $1"_"$2"_"$4"_"$5 }' > ${pop1_folder}/${pop2}/${chr}.${pop1}.tab
                zcat ${inputfile2} | tail -n+2 | awk '{print $1"_"$2"_"$4"_"$5 }' > ${pop2_folder}/${pop1}/${chr}.${pop2}.tab
            ;;
        esac


        #we want the common data
        comm -12 <( sort ${pop1_folder}/${pop2}/${chr}.${pop1}.tab) <( sort ${pop2_folder}/${pop1}/${chr}.${pop2}.tab) > ${outputfolder}/${chr}.${pop1}_${pop2}.common.tab

        case ${format} in
            VCF )
                #now extract this data from each population: we need to generate a new vcf formatted file
                (bcftools view -h ${inputfile} ;bcftools view -H ${inputfile} | awk 'FNR==NR {a[$1]=$0;next} $1"_"$2"_"$4"_"$5 in a {print $0}' ${outputfolder}/${chr}.${pop1}_${pop2}.common.tab - )| bgzip -c > ${common_out_path_pop1}/${chr}.${pop1}_${pop2}.common.vcf.gz
                (bcftools view -h ${inputfile2} ;bcftools view -H ${inputfile2} | awk 'FNR==NR {a[$1]=$0;next} $1"_"$2"_"$4"_"$5 in a {print $0}' ${outputfolder}/${chr}.${pop1}_${pop2}.common.tab - )| bgzip -c > ${common_out_path_pop2}/${chr}.${pop1}_${pop2}.common.vcf.gz
            ;;
            TAB )
                #now extract this data from each population
                zcat ${inputfile} | awk 'FNR==NR {a[$1]=$0;next} $1"_"$2"_"$4"_"$5 in a {print $0}' ${outputfolder}/${chr}.${pop1}_${pop2}.common.tab - | gzip -c > ${common_out_path_pop1}/${chr}.${pop1}_${pop2}.common.tab.gz
                zcat ${inputfile2} | awk 'FNR==NR {a[$1]=$0;next} $1"_"$2"_"$4"_"$5 in a {print $0}' ${outputfolder}/${chr}.${pop1}_${pop2}.common.tab - | gzip -c > ${common_out_path_pop2}/${chr}.${pop1}_${pop2}.common.tab.gz
            ;;
        esac

    ;;
    GENEEXTR)
        # mode=$1
        # inputfile=$2
        # outputfolder=$3
        gene=$4
        cadd_bin=$5

        out_name=`basename ${inputfile}`

        zcat ${inputfile} | awk -v gn=${gene} -v caddbin=${cadd_bin} '{if($8~gn && $(NF-2)==caddbin && $(NF-1) > 0) print $1,$2,$3,$4,$5,gn,caddbin}' > ${outputfolder}/${out_name}.common_enrichment.positive.${cadd_bin}.${gene}.tab
        zcat ${inputfile} | awk -v gn=${gene} -v caddbin=${cadd_bin} '{if($8~gn && $(NF-2)==caddbin && $(NF-1) < 0) print $1,$2,$3,$4,$5,gn,caddbin}' > ${outputfolder}/${out_name}.common_enrichment.negative.${cadd_bin}.${gene}.tab
    ;;
    LOFVCFEXTRLAX )
    # mode=$1
    # inputfile=$2
    # outputfolder=$3

    chr=$4
    pop=$5
    ref_fasta=$6
    # file_type=$6
    # inputfile=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/TSI_chr22_multiSPLIT_LOF.vcf.gz
    # outputfolder=/home/cocca/analyses/paperONE/supp_data/pop_subset/TSI/LOF/22/freq
    input_name=`basename ${inputfile}`

    #define lof variants
    declare -a lof_group
    lof_group=(frameshift_variant splice_acceptor_variant splice_donor_variant start_lost stop_gained stop_lost)
    mkdir -p ${outputfolder}/lof_cons

    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -f ${ref_fasta} -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    # bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    bcftools annotate -x"INFO/AA,INFO/EAS_AF,INFO/AMR_AF,INFO/AFR_AF,INFO/SAS_AF,INFO/EUR_AF" ${inputfile} | bcftools norm -m - | bcftools view -G | bcftools query -f"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC\t%INFO/CADD_PHRED\t%INFO/CSQ\n" | gzip -c > ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz

    # for lof in ${lof_group[@]}
    for (( il = 0; il < ${#lof_group[@]}; il++ ))
    do
        lof=${lof_group[il]}
        # zcat ${inputfile} | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
        # zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
        zcat ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz | awk -v conseq=${lof} '$8~conseq && $8!~"synonymous_variant" && $8!~"missense_variant"' | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
        # zcat ${outputfolder}/${input_name}.maf_csq_bin.tab.gz | fgrep $lof | gzip -c >  ${outputfolder}/lof_cons/${input_name}.maf_${lof}_bin.tab.gz
    done
    #we need also to create a merged version of the data extracted
    zcat ${outputfolder}/lof_cons/${input_name}.maf_*_bin.tab.gz | sort -k1,1 -k2,2 | uniq | gzip -c > ${outputfolder}/lof_cons/${input_name}.maf_all_lof_bin_resume.tab.gz

    #clean splitted file
    rm ${outputfolder}/lof_cons/${input_name}.SPLIT.vcf.gz
    ;;
esac
