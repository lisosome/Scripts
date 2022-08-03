#!/usr/bin/env python3

import argparse
import yaml
import subprocess as sub
import os

#with open("test_config.yaml", "r", encoding="utf8") as f:
    #tst = yaml.safe_load(f)

parser = argparse.ArgumentParser(description='Portable python script to perform CHARGE Depression-Blood pressure analyses')
parser.add_argument('-i', action='store', dest="input_config", help="Input config file in yaml format")
parser.add_argument('-o', action='store', dest="output_path", help="Output path for analyses files")
parser.add_argument('-t', action='store', dest="bp", help="Blood Pressure trait" )
parser.add_argument('-c', action='store', dest='cromo', help='Cromosome for job name')
parser.add_argument('-g', action='store', dest='geno', help='Genotype binary file for MMAP')
parser.add_argument('-s', action='store', dest='sx', help='Sex group to be specified')


params = parser.parse_args()
config_file=params.input_config
gen_file=params.geno
BP=params.bp
out=params.output_path
chr=str(params.cromo)

print("Checking config file existance")
if not os.path.exists(config_file):
    print("Config file doesn't exists. Please provide a config file or check the absolute path")
    exit(1)
else:
    print("Loading config file")
    with open('%s' %(config_file), 'r') as stream:
        config=yaml.safe_load(stream)

print("Checking for output directory...")
if not os.path.isdir(out):
    print("Output directory doesn't exists. Creating a new one")
    os.makedirs(out)
else:
    print("Output directory exists. Proceeding with the analyses...")


sex=params.sx

if config['model']['name'] == "MODEL1":
    if sex == "M":
        if config['model']['depr_trait'] == "dDEPR":
            covs=config['model']['COVs']['dm']
            cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --gxe_interaction {depr} --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{depr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['m'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['m'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
            os.chdir(out)
            sub.check_output(cmd, shell=True, encoding='utf8')
        else:
            covs=config['model']['COVs']['qm']
            cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --gxe_interaction {depr} --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{depr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q fq {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['m'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['m'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
            os.chdir(out)
            sub.check_output(cmd, shell=True, encoding='utf8')

    elif sex == "F":
        if config['model']['depr_trait'] == "dDEPR":
            covs=config['model']['COVs']['df']
            cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --gxe_interaction {depr} --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{depr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q fq {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['f'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['f'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
            os.chdir(out)
            sub.check_output(cmd, shell=True, encoding='utf8')
        else:
            covs=config['model']['COVs']['qf']
            cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --gxe_interaction {depr} --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{depr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q fq {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['f'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['f'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
            os.chdir(out)
            sub.check_output(cmd, shell=True, encoding='utf8')

    elif sex == "comb":
        if config['model']['depr_trait'] == "dDEPR":
            covs=config['model']['COVs']['dcomb']
            cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --gxe_interaction {depr} --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{depr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q fq {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['comb'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['comb'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
            os.chdir(out)
            sub.check_output(cmd, shell=True, encoding='utf8')
        else:
            covs=config['model']['COVs']['qcomb']
            cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --gxe_interaction {depr} --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{depr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q fq {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['comb'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['comb'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
            os.chdir(out)
            sub.check_output(cmd, shell=True, encoding='utf8')

elif config['model']['name'] == "MODEL2":
    if sex == "M":
        covs=config['model']['COVs']['m2']
        cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['m'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['m'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
        os.chdir(out)
        sub.check_output(cmd, shell=True, encoding='utf8')

    elif sex == "F":
        covs=config['model']['COVs']['f2']
        cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['f'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['f'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
        os.chdir(out)
        sub.check_output(cmd, shell=True, encoding='utf8')

    elif sex == "comb":
        covs=config['model']['COVs']['comb2']
        cmd='echo "mmap --ped {ped} --phenotype_filename {pheno} --phenotype_id {id} --trait {press} --covariate_filename {f_cov} --covariates {cov} --read_binary_covariance_file {kin} --single_pedigree --linear_regression --binary_genotype_filename {gen} --hc0_sandwich_estimator --file_suffix {suff} --num_mkl_threads {thr} && gzip *{chr}*.csv" | qsub -N chr{chr}_{press} -j y -o /dev/null -V -l h_vmem={mem} -cwd -q {cal} -pe smp {thr}'.format(ped=config['model']['ped_file'], pheno=config['model']['in_file']['comb'], id=config['model']['pheno_id'], press=BP, mem=config['resources']['mem'], f_cov=config['model']['cov_file']['comb'], cov=covs, kin=config['model']['kinship_file'], depr=config['model']['depr_trait'], gen=gen_file, suff=config['model']['file_suff'] + chr, chr=chr, out=out, thr=config['resources']['threads'], cal=config['resources']['nodes'])
        os.chdir(out)
        sub.check_output(cmd, shell=True, encoding='utf8')
