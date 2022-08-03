#!/bin/bash -f 

# script to convert VCF to MMAP format
# The script uses gunzip to stream the VCF into MMAP function call
# The standard gunzip can be used but gunzip from libdeflate is recommended
# which is 2x faster. A gunzip.libdeflate has been been included, but if
# doesn't run on your machine the code and installation instructions can be found at
# https://github.com/ebiggers/libdeflate
chr=$1
# path to gunzip
gunzip=/home/nardone/software/mmap/mmap_gli_training/programs/gunzip.libdeflate

mmap=/home/nardone/software/mmap/mmap_gli_training/programs/mmap.2022_02_08_14_53.intel

vcf_genotype_file=/netapp06/imputation/IGRPv1/FVG/06022018/03.IMPUTED/VCF/R2/chr${chr}.vcf.gz

#set marker limit to -1 to include all markers
marker_limit=10000
marker_limit=-1

# imputation Rsquared threshold: exlude if < r2 
# recommend at least >0.05
r2=0.00
r2=0.4

# min MAC threshold:  exclude if < total dosage threshold
# if you prefer MAF then solve for corresponding MAC
# The threshold will eliminate very rare variants and speed up conversion 
# recommend 5
total_dosage_threshold=0.0
total_dosage_threshold=5.0

# hard call filter for genotype probabilities
# if the probability > probt then set to 1.0
# recommend 0.99
probt=0.99  

# blocks of markers to stream in bulk
# can be any value >=1
# recommend 1001 
block_size=1001

# mum bits used for encoding genotype probabilities/dosages
# accurate to the 3 decimal place values in VCF  
num_bits=12

# set prefix to output file
mmap_output_prefix="INGI-FVG.vcf_converted.${marker_limit#-}.r2.${r2}.mac.${total_dosage_threshold}.bit.${num_bits}.probt.${probt}.b.${block_size}.chr${chr}"


# ------------------ DO NOT EDIT PAST LINE 

stream_block_size=" --stream_block_size $block_size "
r2_str="--include_r2_field "  # include r2 in converted binary
imputation_str=" "   # do no include IMPUTED/TYPED in binary
imputation_str="--include_imputed_field  "   # include if IMPUTED/TYPED in binary
zstd_compression_level=7  # compression level for ZSTD - 7 generally recommended
num_dosage_bits="--num_dosage_bits ${num_bits} "  # number of bits used to store decimal  
marker_limit_str=" --marker_loop_limit $marker_limit " # limit conversion to number of markers 
imputation_threshold=" --imputation_r2_threshold ${r2} "  #
total_dosage_threshold=" --total_dosage_threshold ${total_dosage_threshold} "
probability_thresholds=" ${probt} ${probt} ${probt} "  # thresholds to set probability to 1.0 (genotype call)
model_str=" --store_GP_genotype_probability  --genotype_probability_best_call_thresholds ${probability_thresholds} "
vcf_commands=" --convert_vcf2mmap  --vcfgz_input_filename ${vcf_genotype_file} "

${gunzip} -cd ${vcf_genotype_file} | \
${dbg} $mmap \
$stream_block_size \
--stream_input \
--mmap_output_prefix ${mmap_output_prefix} \
${min_maf_str} \
${r2_str} \
${imputation_str} \
${vcf_commands} \
${num_dosage_bits} \
${imputation_threshold}  \
${total_dosage_threshold}  \
--use_zstd \
--zstd_compression_level ${zstd_compression_level}  \
${marker_limit_str}  \
${model_str}
