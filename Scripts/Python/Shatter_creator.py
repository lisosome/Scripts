#!/usr/bin/env python3
#from pysam import VariantFile as vf
import sys
import os
import pandas as pd
import re
import gzip

infile=sys.argv[1]
outpath=sys.argv[2]
#vcfin=vf(infile)
#dup=re.compile('chr\d{1,2}:\d{1,10}')
dup=re.compile('MantaDUP')
ins=re.compile('MantaINS')
dl=re.compile('MantaDEL')
bnd=re.compile('MantaBND')
#for rec in vcfin.fetch():
    #if dup.search(rec.id):
        #with vf(outpath +"/"+"Tandem_dup.txt", "w") as out:
            #out.write(rec)
            #out.close()
    #elif ins.search(rec.id):
        #with vf(outpath +"/"+"Insertions.txt", "w") as inser:
            #inser.write(rec)
            #inser.close()
    #elif dl.search(rec.id):
        #with vf(outpath +"/"+"Dels.txt", "w") as dl:
            #dl.write(rec)
            #dl.close()
    #elif bnd.search(rec.id):
        #with vf(outpath +"/"+"Tandem_dup", "w") as bd:
            #bd.write(rec)
            #bd.close()
for line in gzip.open(infile, 'r'):
    line=line.decode('utf8')
    if not(re.match('#', line.strip())):
        id=line.strip().split("\t")[2]
        info=line.strip().split("\t")[7]
        #head="chrom1\tpos1\tchrom2\tpos2\tSVtype\tstrand1\tstrand2"
        if dup.search(id):
            chr=line.strip().split("\t")[0]
            ps1=int(line.strip().split("\t")[1])
            ps2=ps1+1
            end1=int(info.split(";")[0].split("=")[1])
            end2=end1+1
            strand1="-"
            strand2="+"
            typ="DUP"
            #var_frame=pd.DataFrame({'chrom1':[],'pos1':[],'chrom2':[],'pos2':[],'SVtype':[],'strand1':[],'strand2':[]})
            #var_frame=pd.DataFrame(pd.concat([chr,ps1,chr,ps2,typ,strand1,strand2]))
            with open(outpath +"/"+"Tandem_dup.txt", "a+") as out:
                #pd.set_option('max_colwidth', None)
                #out.write(head+"\n")
                out.write(chr+"\t"+str(ps1)+"\t"+chr+"\t"+str(end2)+"\t"+typ+"\t"+strand1+"\t"+strand2+"\n")
                out.close()
        #elif ins.search(id):
            #with open(outpath +"/"+"Insert.txt", "a+") as out:
                #out.write(line)
                #out.close()
        #elif dl.search(id):
            #with open(outpath +"/"+"Del.txt", "a+") as out:
                #out.write(line)
                #out.close()
        elif bnd.search(id):
            with open(outpath +"/"+"Breakend.txt", "a+") as out:

                out.write(line)
                out.close()
