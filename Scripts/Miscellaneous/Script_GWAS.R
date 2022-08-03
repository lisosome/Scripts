## calcolo kinship con GEMMA
/shared/software/gemma/gemma-0.98-linux-static  -bfile /netapp02/data/genotypes/SR/SR_updated_20210312 -gk 2 -o /home/nardone/SR_updated_20210312_kinship
##

##regressione con GEMMA
chr=22
geno=/netapp02/data/imputation/HRC/SR/BIMBAM/chr${chr}.pbwt_reference_impute_clean.gen.gz.bimbam.gz
pheno=/home/nardone/File_GEMMA/PHENO_893_D_P_SR.txt
snp_ann=/netapp02/data/imputation/HRC/SR/BIMBAM/chr${chr}.pbwt_reference_impute_clean.gen.gz.pos
kinship=/home/nardone/output/SR_updated_20210312_kinship.sXX.txt
cov=/home/nardone/File_GEMMA/Covariate_sesso_eta_scolarita_D_P_SR.txt
out=GEMMA_SR_discromie_D_P_chr${chr}.txt
echo "/shared/software/gemma/gemma-0.98.4-linux-static-AMD64-g ${geno} -p ${pheno} -a ${snp_ann} -k ${kinship} -lm 4 -maf 0 -c ${cov} -o ${out}" |  qsub -N job_${chr}_v  -o \$JOB_ID_${chr}_v.log -e \$JOB_ID_${chr}_v.e -V -l h_vmem=10G -cwd -q fast,all.q@geneticalabdb7.burlo.trieste.it

## for loop regressione Gemma
for chr in {1..21}
do
geno=/netapp02/data/imputation/HRC/SR/BIMBAM/chr${chr}.pbwt_reference_impute_clean.gen.gz.bimbam.gz
pheno=/home/nardone/File_GEMMA/PHENO_893_D_P_SR.txt
snp_ann=/netapp02/data/imputation/HRC/SR/BIMBAM/chr${chr}.pbwt_reference_impute_clean.gen.gz.pos
kinship=/home/nardone/output/SR_updated_20210312_kinship.sXX.txt
cov=/home/nardone/File_GEMMA/Covariate_sesso_eta_scolarita_D_P_SR.txt
out=GEMMA_SR_discromie_D_P_chr${chr}.txt
echo "/shared/software/gemma/gemma-0.98-linux-static -g ${geno} -p ${pheno} -a ${snp_ann} -k ${kinship} -lm 4 -maf 0 -c ${cov} -o ${out}" |  qsub -N job_${chr}_v  -o \$JOB_ID_${chr}_v.log -e \$JOB_ID_${chr}_v.e -V -l h_vmem=10G -cwd -q fast,all.q@geneticalabdb7.burlo.trieste.it
done

##selezione 0.01<MAF<0.99 [Ubuntu]
chr=22
echo -e 'chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score' > /home/nardone/Discromie/Deutano_Protano/File_Filt_MAF01/SR_D_P_maf01_chr${chr}.txt
cat /home/nardone/Discromie/Deutano_Protano/output/GEMMA_SR_discromie_D_P_chr22.txt.assoc.txt | sed '1d' |  awk 'BEGIN {FS=" "} ($8>=0.01 && $8<=0.99) {print $0}' >> /home/nardone/Discromie/Deutano_Protano/Test_MAF/SR_D_P_maf01_chr${chr}.csv

## for loop MAF [Ubuntu]
for chr in {1..21}
do
echo -e 'chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score'> /home/nardone/Discromie/Deutano_Protano/File_Filt_MAF01/SR_D_P_maf01_chr${chr}.txt
cat /home/nardone/Discromie/Deutano_Protano/output/GEMMA_SR_discromie_D_P_chr${chr}.txt.assoc.txt | sed '1d' |  awk 'BEGIN{FS=" "} ($8>=0.01 && $8<=0.99) {print $0}' >> /home/nardone/Discromie/Deutano_Protano/File_Filt_MAF01/SR_D_P_maf01_chr${chr}.txt
done

## recupero RSQ per Silk Road
cd /netapp02/data/imputation/HRC/SR/prob
file .info

######## in R #####
#codice per recuperare rsq in R
chr=22
file_name<-paste("SR_D_P_maf01_chr",chr,".txt",sep="")
d<-read.table(file=file_name,he=T)
dim(d) # 5678879      13
info<-read.table(file=paste("/netapp02/data/imputation/HRC/SR/prob/SR_HRC_chr",chr,".info",sep=""),he=T)
info<-info[info$Rsq>=0.4,]
## Aggiungere colonna rs_allele1_allele0 nel file info per poter confrontare i file
info$key<-paste(info$rs_id,info$A0,info$A1,sep="_")
##Unire i due file
d1<-d[d$rs%in%info$key,] #prendo righe in comune tra i due files
write.table(d1,"/home/nardone/Discromie/Deutano_Protano/Test_RSQ/SR_Discromie_rsid_rsq04_maf01.txt",quote=F,col.names = T,row.names = F, sep="\t")

##Loop automatizzazione step precedenti creando file che comprendono gli rs in R
for(chr in 1:21){
	file_name<-paste("SR_D_P_maf01_chr",chr,".txt",sep="")
	d<-read.table(file=file_name,he=T)
	info<-read.table(file=paste("/netapp02/data/imputation/HRC/SR/prob/SR_HRC_chr",chr,".info",sep=""),he=T)
	 info<-info[info$Rsq>=0.4,]
	 info$key<-paste(info$SNP,info$A0,info$A1,sep="_")
	 d1<-d[d$rs%in%info$key,]
	 file_out=paste("/home/nardone/Discromie/Deutano_Protano/Test_RSQ/chr", chr,"_SR_Discromie_rsid_rsq04_maf01.txt", sep="")
	 write.table(d1,file_out,col.names = T,row.names = F, sep="\t")
}
## Loop per senza rs in R
for(chr in 1:22){
	file_name<-paste("SR_D_P_maf01_chr",chr,".txt",sep="")
	d<-read.table(file=file_name,he=T)
	info<-read.table(file=paste("/netapp02/data/imputation/HRC/SR/prob/SR_HRC_chr",chr,".info",sep=""),he=T)
	 info<-info[info$Rsq>=0.4,]
	 d1<-d[d$rs%in%info$SNP,]
	 file_out=paste("/home/nardone/Discromie/Deutano_Protano/Test_RSQ/chr", chr,"_SR_Discromie_NOrs_rsq04_maf01.txt", sep="")
	 write.table(d1,file_out,col.names = T,row.names = F, sep="\t")
}

### da linea di comando ####
#uniamo i files

echo -e 'chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score' > /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR.txt
for chr in {1..22}
do
cat /home/nardone/Discromie/Deutano_Protano/File_MAF01_RSQ04/chr${chr}_SR_Discromie_rsid_rsq04_maf01.txt | sed '1d' | awk 'BEGIN {FS=" "}  {print $0}' >> /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR.txt
cat /home/nardone/Discromie/Deutano_Protano/File_MAF01_RSQ04/chr${chr}_SR_Discromie_NOrs_rsq04_maf01.txt | sed '1d' | awk 'BEGIN {FS=" "}  {print $0}' >> /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR.txt
done

##Fare file best selezionando per il p-value [Ubuntu]
echo -e 'chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score' > /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR_BEST.txt
cat /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR.txt | sed '1d' | awk 'BEGIN {FS=" "}  ($13<=0.0001) {print $0}' >> /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR_BEST.txt

## VEP da riga di comando
/shared/software/VEP/ensembl-vep-100/ensembl-vep/vep  -i /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/FILE_BEST_DISCRO_perVEP.txt --cache --dir_cache /shared/software/VEP/ensembl-vep/modules/Bio/EnsEMBL/VEP/Config.pm -o /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/VEP_BEST_DISCRO.txt -e --distance [250000]
###############################################

##Manhattan e QQ plot

screen -S grafici_GWAS  #creo screen separato in cui allochiamo 10 giga di memoria
qrsh -V -l h_vmem=10G

#Diminuire quantità di dati: dal file ALL_CHROMOSOME tolgo un po' di SNP per ridurre le dimensioni dei plot
data2<-a[a$p_score>=0.05,]
  x<-sample(rownames(a),500000)
  data_ns<-a[rownames(a)%in%x,]
  data_s<-a[a$p_score<0.05,]
  data_ok<-rbind(data_ns,data_s)



#Script Manhattan
pdf("Manhattan_SR_D_P_.pdf",width=21)
  manhattan(data_ok, chr = "chr", bp = "ps", p = "p_score", snp = "rs", col = c("orange", "blue"), chrlabs = NULL, suggestiveline = -log10(1e-06), genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE,cex.main=2.5)
  dev.off()
##Script QQ plot
pdf("QQplot_SR_D_P_MAF05.pdf")
qq(a$p_score)
dev.off()
## Calcolo lambda
lambda <- round(median((a$beta/a$se)^2)/ 0.454, 3)
lambda

echo -e 'chr\trs\tps\tn_mis\tn_obs\tallele1\tallele0\taf\tbeta\tse\tp_wald\tp_lrt\tp_score' > /home/nardone/Discromie/GWAS_Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR_BEST_MAF05.txt
cat /home/nardone/Discromie/GWAS_Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR.txt | sed '1d' | awk 'BEGIN {FS=" "}  ($8>=0.05 && $8<=0.95) {print $0}' >> /home/nardone/Discromie/Deutano_Protano/File_ALL_CHROMOSOME/SR_D_P_ALLCHR_BEST_MAF05.txt


##Principal Componnent Analysis (PCA) con plink
plink --bfile /netapp02/data/genotypes/SR/SR_updated_20210312 --indep-pairwise 50 5 0.5 ##filtrare files per r^2=0.5
plink --bfile /netapp02/data/genotypes/SR/SR_updated_20210312 --extract plink.prune.in --make-bed --out Prunedata_SR ##fare file con SNP filtrati
plink --bfile Prunedata_SR --pca 10 --out PCA_SR ##comando PCA
