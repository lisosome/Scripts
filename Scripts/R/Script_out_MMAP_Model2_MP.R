args <- commandArgs(trailing=TRUE)
chr <- args[[1]]

############ ****** QC e preparazione file output interaction ****** ######### 
library(stringr)
#chr=22
STUDY="INGI-FVG"
IMP="IGRPv1"
PHENO="QT"
MODEL="M2"
DATA="20190515"

x<-paste("apro file ",PHENO,".model2_",chr,".lm.pval.csv.gz",sep="")
write.table(x,paste("log_mio_chr",chr,".txt",sep=""))
d<-read.csv(file=paste(PHENO,".model2_",chr,".lm.pval.csv",sep=""),he=T)

x<-"elimino le colonne che non servono"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
d1<-d[,c(1:6,9,11,13,17,18,21)]
rm(d)

# faccio un primo QC per maf e NA
x<-"QC per maf e NA in file output"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
d1<-d1[!is.na(d1$POS),]
d1<-d1[which(d1$EFFECT_ALLELE_FREQ>=0.01 & d1$EFFECT_ALLELE_FREQ<=0.99),]
d1<-d1[which(!is.na(d1$SNP_TSCORE_PVAL)),]

## apro info
x<-"apro info"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

info<-read.table(file=paste("/netapp02/data/imputation/INGI_TGP3_MERGED/CARL/MERGED/ALL/RELEASE/CONVERTED/chr",chr,".gen_info",sep=""),he=T)
# QC per info e maf
x<-"QC per info e maf in file info"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

info<-info[which(info$info>=0.1),]
info<-info[which(info$exp_freq_a1>=0.01 & info$exp_freq_a1<=0.99),]

x<-"creo colonna segn in entrambi i file"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

info$segn<-paste(info$rs_id,info$a0,info$a1,sep="_")
d1$segn<-paste(d1$RSNUM,d1$NON_CODED_ALLELE,d1$EFFECT_ALLELE,sep="_")

x<-"controllo intersezione info e output"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
d2<-d1[which(d1$segn%in%info$segn),]  
info2<-info[which(info$segn%in%d1$segn),]
info3<-info[which(!(info$segn%in%d1$segn)),]
d3<-d1[which(!(d1$segn%in%info$segn)),] 

x<-"segn per ins/del in info"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

info3$segn<-paste(info3$rs_id,info3$a1,info3$a0,sep="_")
info4<-info3[which(info3$segn%in%d3$segn),]
info5<-info3[which(!(info3$segn%in%d3$segn)),]
info_ok<-rbind(info2,info4)


info5$segn<-paste(info5$rs_id,"R_D",sep="_")
info4<-info5[which(info5$segn%in%d1$segn),]
info5<-info5[which(!(info5$segn%in%d1$segn)),]
info_ok<-rbind(info_ok,info4)

info5$segn<-paste(info5$rs_id,"R_I",sep="_")
info4<-info5[which(info5$segn%in%d1$segn),]
info5<-info5[which(!(info5$segn%in%d1$segn)),]
info_ok<-rbind(info_ok,info4)


x<-"aggiorno output file e controllo per valori multipli d segn"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

d3<-d1[which(d1$segn%in%info_ok$segn),]
d2<-d1[which(!(d1$segn%in%info_ok$segn)),] ## sono quelli con info inferiore a 0.1
a<-table(d3$segn)
b<-a[a>1]
d4<-d3[which(!(d3$segn%in%names(b))),]

x<-"aggiorno info file e tolgo multipli di segn"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

info_ok<-info_ok[which(info_ok$segn%in%d4$segn),]
a<-table(info_ok$segn)
b<-a[a>1]
info_ok1<-info_ok[which(!(info_ok$segn%in%names(b))),]

x<-"aggiorno output file"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

d5<-d4[which(d4$segn%in%info_ok1$segn),]

x<-"controllo che i due file abbiano lo stesso numero di marker"
write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)

if(dim(d5)[1]==dim(info_ok1)[1]){
  x<-"info e output stesso numero di snp... procedo con il merge"
  write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
  OK<-merge(d5,info_ok1,by="segn")
  x<-paste("numero totale di SNP e colonne ",dim(OK)[1]," ",dim(OK)[2],sep="")
  write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
  x<-"scrivo file di output finale"
  write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
  OK1<-OK[,c(2:13,20,22)]
  
  write.table(OK1,paste(STUDY,".EA.",IMP,".",PHENO,".CALCIUM",".",MODEL,".chr",chr,".txt",sep=""),quote=F,col.names = T,row.names = F)
  
  
  x<-summary(OK1)
  write.table(x,paste("summary_mio_chr",chr,".txt",sep=""))
  
  x<-"FINITO ----> TUTTO OK!!!!"
  write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
  
}
if(dim(d5)[1]!=dim(info_ok1)[1]){
  x<-"info e output numero diverso di snp... ATTENZIONE --> NIENTE MERGE"
  write.table(x,paste("log_mio_chr",chr,".txt",sep=""),append=T)
}


### lancio su bash con code
#for chr in {1..20}
#do
 
#echo "R CMD BATCH '--args '${chr}' ' /netapp04/concas/dose_MMAP/FVG/interaction/script_sistemazioneOutputMMAP.R " | qsub -N job${chr}  -o \$JOB_ID_${chr}.log -e \$JOB_ID_${chr}.e  -V -l h_vmem=30G -cwd -q fast
 
#done      



### file unico 
# echo 'SNPNAME RSNUM CHR POS NON_CODED_ALLELE EFFECT_ALLELE NUM_OBS EFFECT_ALLELE_FREQ BETA_SNP SNP_TSCORE_PVAL HC0_SE_SNP SNP_CHI_SQR_PVAL BETA_G_x_CALCIUM PVAL_G_x_CALCIUM HC0_SE_G_x_CALCIUM HC0_COV_BETA_CALCIUM_BETA_G HC0_COV_BETA_CALCIUM_BETA_G_x_CALCIUM HC0_COV_BETA_G_BETA_G_x_CALCIUM info type' > INGI-CAR.EA.IGRPv1.QT.CALCIUM.M1.20190510.txt
# 
# for chr in {1..22}
# do
# cat INGI-CAR.EA.IGRPv1.QT.CALCIUM.M1.chr${chr}.txt | sed '1d' | awk 'BEGIN {FS=" "}  {print $0}' >> INGI-CAR.EA.IGRPv1.QT.CALCIUM.M1.20190510.txt
# done

