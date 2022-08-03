#!/usr/bin/env bash

#Let's storm some ideas to execute SV variant calling with manta and delly. Then we will merge with survivor
#In order to make manta run, a running script must be generated from the configuration one. But first, let's define some folders
#Since the crams output will be provided from Max's pipeline, we have to create another manifest. (this can be avoided if we incorporate these procedures in Max's ppipeline, but that's a story for another time)


#Define a master output in which all the subfolders of the samples will be created
Mout=$1 #master output folder. Something like 3.SV_VARCALL
#Rather than provide an alredy created manifest, why don't wa make the script do it! Yay! We even don't need a header. But we need a sample list file (bummer)
Slist=$2 #list of samples. One sample per line
Fpath=$3 #path in which are located bam/cram to analyze. From Max's pipeline it will be something like /large/___SCRATCH___/burlo/user/folder/VAR_CALL_OUT_FOLDER/2.BQSR

#Let's verify that our beloved and needed conda environment are present
mkdir -p ${Mout}

module load conda samtools bcftools
ref_genome=/fast/burlo/fcrudele/resources/hgRef/GRCh38.p13/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna

while read line
do
  ID=$(echo $line)
  cram=$(find -L ${Fpath} -type f -name "${ID}_bqsr.cram")
  paste -d " " <(echo "${ID}") <(echo "${cram}") >> ${Mout}/SV_manifest.tsv
done < ${Slist}

conda_u=/large/___HOME___/burlo/nardone/.conda/envs
#Let's calculate coverage
coverOut=${Mout}/coverage
bed=${coverOut}/bedFiles
#callable=${coverOut}/callableIntervals
#new_bams=${coverOut}/Converted_BAMs
mkdir -p ${coverOut}
mkdir -p ${bed}
#let's export some usefule env vars
export MOSDEPTH_Q0=NO_COVERAGE   # -- defined by the arguments to --quantize
export MOSDEPTH_Q1=LOW_COVERAGE
export MOSDEPTH_Q2=CALLABLE
export MOSDEPTH_Q3=HIGH_COVERAGE

conda activate smoove_pers
#Like this, we generate a gzipped bed file per sample, binned for coverage values and filter callable intervals. Then, we filter cram files
while read line
do
  sam=$(echo $line | awk '{print $1}')
  cram=$(echo $line | awk '{print $2}')
  out_file=${bed}/${sam}_coverage.bed
  mosdepth -n --quantize 0:1:4:150: -f ${ref_genome} ${out_file} ${cram}
  #sam_fold=${new_bams}/${sam}
  #mkdir -p ${sam_fold} && samtools view -T ${ref_genome} -b -o ${sam_fold}/${sam}_converted.bam ${cram} && samtools index ${sam_fold}/${sam}_converted.bam
done < ${Mout}/SV_manifest.tsv

echo "Coverage calculation completed. Now let's call some structural variants"

#First, let's update our lovable manifest
#while read line
#do
  #sam=$(echo $line | awk '{print $1}')
  #sam_fold=${new_bams}/${sam}
  #old_path=$(echo $line | awk '{print $2}')
  #new_path=$(find -L ${sam_fold} -type f -name "*.bam")
  #paste -d " " <(echo "${sam}") <(echo "${old_path}") <(echo "${new_path}") >> ${Mout}/SV_manifest.tsv
#done < ${Mout}/SV_manifest.tsv
conda deactivate
#This is were the fun begins. Let's call manta environment
#module load conda
conda activate Manta

#Let's define some variables for easy coding
PathMan=${conda_u}/Manta/bin
#We generate Manta's run scripts through the configuration one
while read line
do
sam=$(echo $line | awk '{print $1}')
bam=$(echo $line | awk '{print $2}')
python2.7 ${PathMan}/configManta.py --bam ${bam} --referenceFasta ${ref_genome} --runDir ${Mout}/1.Manta_outs/${sam}/Manta_run
#let's call the SVs!
res=${Mout}/1.Manta_outs/${sam}/Manta_run/results/variants/diploidSV.vcf.gz
samtools=/large/___HOME___/burlo/nardone/.conda/envs/Manta/share/manta-1.6.0-1/libexec/samtools
python2.7 ${Mout}/1.Manta_outs/${sam}/Manta_run/runWorkflow.py
python2.7 /large/___HOME___/burlo/nardone/.conda/envs/Manta/share/manta-1.6.0-1/libexec/convertInversion.py ${samtools} ${ref_genome} ${res} > ${Mout}/1.Manta_outs/${sam}/${sam}_final.vcf
done < ${Mout}/SV_manifest.tsv
conda deactivate

echo "Manta's calling completed. Proceeding with delly"

conda activate del9

delly=${conda_u}/del9/delly_v0.9.1_linux_x86_64bit
#map=/fast/burlo/fcrudele/resources/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz
delout=${Mout}/2.Delly_outs
mkdir -p ${delout}
while read line
do
  sam=$(echo $line | awk '{print $1}')
  bam=$(echo $line | awk '{print $2}')
  ${delly} call -g ${ref_genome} -o ${delout}/${sam}.bcf ${bam}
done < ${Mout}/SV_manifest.tsv

#Now let's define a variable with all the bcf samples names
bcf=$(echo ${delout}/*.bcf)
#merging step
${delly} merge -o ${delout}/merged_sites.bcf ${bcf}
#Let's genotype
while read line
do
  sam=$(echo $line | awk '{print $1}')
  bam=$(echo $line | awk '{print $2}')
  ${delly} call -g ${ref_genome} -v ${delout}/merged_sites.bcf -o ${delout}/${sam}_geno.bcf ${bam}
done < ${Mout}/SV_manifest.tsv

conda deactivate
#One last effort. Let's convert all the bcfs in vcfs. But let's do it in a separate folder
vcfs=${delout}/delly_VCFs
filter=${vcfs}/filtered
mkdir -p ${vcfs}
mkdir -p ${filter}
while read line
do
  sam=$(echo $line | awk '{print $1}')
  bcftools view ${delout}/${sam}_geno.bcf > ${vcfs}/${sam}_geno.vcf &&
  bcftools view -i "FILTER=='PASS'" ${vcfs}/${sam}_geno.vcf -Ov -o ${filter}/${sam}_filtered.vcf
done < ${Mout}/SV_manifest.tsv

echo "Delly calling completed. Almost there. Now let's Smoove things up!"

conda activate smoove_pers

skip_bed=/large/___SCRATCH___/burlo/nardone/Smoove_test/test_out_module/test_samples/regions_to_skip.bed
lenght=$(wc -l ${Mout}/SV_manifest.tsv )
readarray -t bams < <(awk '{print $2}' ${Mout}/SV_manifest.tsv)
bam_p=$(echo "${bams[@]}")
sm_outs=${Mout}/3.Smoove_outs
smVCF=${sm_outs}/Smoove_VCFs
mkdir -p ${sm_outs}
mkdir -p ${smVCF}
if [[ ${lenght} < 40 ]]
then
smoove call -x --name smoove_calling --outdir ${sm_outs} --exclude ${skip_bed} --fasta ${ref_genome} -p 3 --genotype ${bam_p}
#Let's split in single sample VCFs
while read line
do
  sam=$(echo $line | awk '{print $1}')
  bcftools view -s ${sam} ${sm_outs}/smoove_calling-smoove.genotyped.vcf.gz -Ov -o ${smVCF}/${sam}_smoove.vcf &&
  rm ${sm_outs}/${sam}*.bam
done < ${Mout}/SV_manifest.tsv
else
  echo "WARNING: analyzed cohort contains more than ~40 samples. This script currently doesn't support all smoove workflow"
  conda deactivate
  exit 1
fi


conda deactivate

echo "We survived so far. So let's use SURVIVOR to merge the outputs!"
conda activate survivor

#first thing first, we need to create n files with the paths of the SV_called files to merge.
merlis=${Mout}/merging_lists
mkdir -p ${merlis}
while read line
do
  sam=$(echo $line | awk '{print $1}')
  find -L ${Mout}/1.Manta_outs/${sam} -type f -name "*${sam}_final.vcf" >> ${merlis}/${sam}_files_to_merge.txt
  find -L ${filter} -type f -name "${sam}_filtered.vcf" >> ${merlis}/${sam}_files_to_merge.txt
  find -L ${smVCF} -type f -name "${sam}_smoove.vcf" >> ${merlis}/${sam}_files_to_merge.txt
done < ${Mout}/SV_manifest.tsv

#Then, we need to merge the specified files using survivor
sur=${Mout}/4.Survivor/merge
stat=${Mout}/4.Survivor/stats
#plots=${Mout}/4.Survivor/plots
mkdir -p ${sur} ${stat}
while read line
do
  sam=$(echo $line | awk '{print $1}')
  SURVIVOR merge ${merlis}/${sam}_files_to_merge.txt 1000 2 1 1 0 30 ${sur}/${sam}_merged.vcf &&
  SURVIVOR stats ${sur}/${sam}_merged.vcf -1 -1 -1 ${stat}/${sam}_stats
  bcftools sort ${sur}/${sam}_merged.vcf -Ov -o ${sur}/${sam}_unnotated.vcf &&
  rm ${sur}/${sam}_merged.vcf
  #SURVIVOR genComp ${sur}/${sam}_merged.vcf ${sam}_merged_matr.txt
done < ${Mout}/SV_manifest.tsv
conda deactivate

echo "Files succesfully merged! Now, let's finish with the annotation"
ann=${Mout}/5.Annotated
mkdir -p ${ann}
cod=${ann}/coding
#rare=${ann}/rare
#rco=${rare}/coding
mkdir -p ${cod}
#mkdir -p ${rco}
svpack=/large/___HOME___/burlo/nardone/svpack/svpack
gff=/large/___HOME___/burlo/nardone/svpack/resources/ensembl.GRCh38.101.reformatted.gff3
#cntr=/large/___HOME___/burlo/nardone/svpack/resources/HPRC_GIAB.GRCh38.pbsv.vcf.gz
#n order to make svpack operate, we have to have the pysam module. Luckily we have it in a conda env
conda activate pysam
#Let's annotate the files with svpack in general, then extracting only variants affecting coding sequences and than compare our results with a set of controls
while read line
do
sam=$(echo $line | awk '{print $1}')
${svpack} consequence ${sur}/${sam}_unnotated.vcf /large/___HOME___/burlo/nardone/svpack/resources/ensembl.GRCh38.101.reformatted.gff3 > ${ann}/${sam}_annotated.vcf
done < ${Mout}/SV_manifest.tsv

while read line
do
sam=$(echo $line | awk '{print $1}')
${svpack} consequence --require-csq ${sur}/${sam}_unnotated.vcf /large/___HOME___/burlo/nardone/svpack/resources/ensembl.GRCh38.101.reformatted.gff3 > ${cod}/${sam}_ann_coding_only.vcf
done < ${Mout}/SV_manifest.tsv

conda deactivate

echo "Workflow completed. Results available at ${ann} and ${cod}"
