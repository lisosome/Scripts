#!/usr/bin/env bash

#PBS -N smoove-test
#PBS -q thin
#PBS -V
#PBS -k eod
#PBS -l select=1:ncpus=5:mem=15g
#PBS -l walltime=96:00:00
#PBS -J 0-16
#PBS -o /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/Log/geno_log/smoove-script_geno_${PBS_ARRAY_INDEX}.log
#PBS -e /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/Log/geno_log/smoove-script_geno_${PBS_ARRAY_INDEX}.err

skip_bed=/large/___SCRATCH___/burlo/nardone/Smoove_test/test_out_module/test_samples/regions_to_skip.bed
reference_fasta=/storage/burlo/cocca/resources/hgRef/hg38/Homo_sapiens_assembly38.fasta
cram_path=/large/___SCRATCH___/burlo/nardone/Smoove_test/cram_samples

module load conda bcftools samtools htslib gnu/9.3.0
conda activate smoove_pers

OPTIONAL_INPUT=${cram_path}/smoove_samples_${PBS_ARRAY_INDEX}.txt

#calling genotype
#while read line
#do
#name=$(echo -e $line)
#smoove call -d --outdir /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/call_results --exclude ${skip_bed} --name ${name} --fasta ${reference_fasta} -p 1 --genotype ${cram_path}/${name}_bqsr.cram &&
#rm /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/call_results/${name}*.bam
#done < ${OPTIONAL_INPUT}

#merging
#smoove merge --outdir /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/merge --name test_821_geno -f ${reference_fasta} /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/call_results/*.genotyped.vcf.gz


#genotyping
while read line
do
sam=$(echo -e $line)
smoove genotype -d -x -p 1 --outdir /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/geno_results --name ${sam}-joint --fasta ${reference_fasta} --vcf /large/___SCRATCH___/burlo/nardone/Smoove_test/pers_env_/test_821/merge_results/821_merge.sites.vcf.gz ${test_path}/${sam}_bqsr.cram
done < ${OPTIONAL_INPUT}
