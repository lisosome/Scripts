#!/usr/bin/env python3

import sys
import os
import subprocess as sub
import dask.dataframe as pd
import gzip
import io
import re


input_file=sys.argv[1]

#sys.stderr=io.StringIO()

if input_file == "-h":
    print("\nUse this script to serach for variants within regions of interest in WGS or imputed data from all cohorts.")
    print("Regions must be contained in a bed file, with information for chromosome, starting position and ending position.\nColumns MUST be tab separated\n\n***ATTENTION***\nChromosome must be specified as 1 or chr1 respectively for GRCh37 and GRCh38")
    print("\nPopulations available for the search:\n\nWGS\n\tGRCh37\n\t\tFVG: Low coverage WGS data of FVG cohort\n\t\tCAR: Low coverage WGS data of Carlantino cohort\n\t\tVBI: Low coverage WGS data of Val Borbera cohort\n\tGRCh38\n\t\tSR: High coverage WGS data from the Silk Road cohort\n\t\tHCFVG: High coverage WGS data from the FVG cohort\n\t\tCOV: High coverage WGS data from the Covid19 cohort\n\t\tOTO: High coverage WGS data from the Otorino cohort\n\t\tWCS: High coverage WGS data from all the HCFVG, COV and OTO callset")
    print("\nImputed data (IMP)\n\tMOL: Imputed data from the Molisani cohort\n\tBAU: Imputed data from the BAU cohort\n\tCAR: Imputed data from the Carlantino cohort\n\tDIABETES: Imputed data from the Diabetes cohort\n\tERB: Imputed data from the Erborinati cohort\n\tI_FVG: Imputed data from the FVG cohort\n\tITT: Imputed data from the Italian Taste cohort\n\tNUTRIACT: Imputed data from the Nutriact cohort\n\tSLO: Imputed data from the Slovenian cohort\n\tVBI: Imputed data from the Val Borbera cohort")
    print("\nSelect between PROMPT and JOB mode. The first in an interactive mode thought to be more user friendly.\nThe latter is more advanced, thought to be used in iterations or to launch the script as a job.\nIn JOB mode, present the needed arguments in the follow order, like a bash script:\n\n1. Bed file\n2. build of reference\n2.5 type of data, only for GRCh37 reference\n3. population\n4. personal name, a name tag to personalize the output file")
    print("\ne.g JOB mode for build 37: ./pop_variant_screen_v1.py example.bed JOB 37 WGS FVG name_to_add\n    JOB mode for build 38: /pop_variant_screen_v1.py example.bed JOB 38 HCFVG name_to_add")
    exit(1)
mode=sys.argv[2]

pop_path37={"FVG":"/shared/INGI_WGS/GRCh37/FVG", "CAR":"/shared/INGI_WGS/GRCh37/CAR", "VBI":"/shared/INGI_WGS/GRCh37/VBI"}
pop_path38={"SR":"/shared/WGS/SilkRoad/20210322_RELEASE_FILTERED", "HCFVG":"/netapp06/WGS/RELEASE/20220210/GRCh38/FVG/1.PHASED", 'COV':"/netapp06/WGS/RELEASE/20220210/GRCh38/COVIDALL/1.PHASED", 'OTO': "/netapp06/WGS/RELEASE/20220210/GRCh38/OTORINO/1.PHASED", 'WCS':"/netapp06/WGS/RELEASE/20220210/GRCh38/WHOLE_CALLSET/1.PHASED", 'OTO2':'/netapp06/WGS/RELEASE/20220505/OTORINO/RECALC', 'HCFVG2':'/netapp06/WGS/RELEASE/20220505/FVG'}
imputed_data={'MOL':'/netapp06/imputation/IGRPv1/MOLISANI/20210613/06.imputed/VCF', 'BAU':'/netapp06/imputation/IGRPv1/BAU/20211210/03.IMPUTED/VCF','CAR':'/netapp06/imputation/IGRPv1/CAR/06022018/03.IMPUTED/VCF/R2','DIABETES':'/netapp06/imputation/IGRPv1/DIABETES/20211210/03.IMPUTED/VCF','ERB':'/netapp06/imputation/IGRPv1/ERBORINATI/20211210/03.IMPUTED/VCF', 'I_FVG':'/netapp06/imputation/IGRPv1/FVG/06022018/03.IMPUTED/VCF/R2','ITT':'/netapp06/imputation/IGRPv1/ITT/20211209/B1_B2_B303.IMPUTED/  VCF', 'NUTRIACT':'/netapp06/imputation/IGRPv1/NUTRIACT/20211210/03.IMPUTED/VCF', 'SLO':'/netapp06/imputation/IGRPv1/SLOVENI/20211210/03.IMPUTED/VCF', 'VBI':'/netapp06/imputation/IGRPv1/VBI/06022018/03.IMPUTED/VCF/R2'}
CHR=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,"X")
#data_type=str(input("Please insert the type of data you want to query.\n\nChoose between WGS or IMP (imputed data): "))
#regions=pd.read_table(input_file, sep="\t", names["cro", "start", "end"])
regions=pd.read_csv(input_file, sep="\t", names=["cro", "start", "end","gene"], dtype='object', engine='c')


if mode == "PROMPT":
    build_input=int(input("Please insert the desired reference build\nChoose between:\n37 or 38: "))

    def bcf_search(vcf):
        if  pop_input != "SR":
            split_fields="%ID\\t%SYMBOL\\t%Feature\\t%Consequence\\t%gnomAD_AF\\t%CDS_position\\t%EXON\\t%INTRON\\t%Protein_position\\t%Amino_acids\\t%SIFT\\t%PolyPhen\\t%CADD_RAW[\\t%GT]\\n"
            #awk_carr="""awk 'BEGIN{FS="\\t"}{for(i=1;i<14;i++) {printf "%s ",$i };for(j=14;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",$j};printf "\\n";}' """
            #tr='tr " " "\t"'
            query='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\t%MAF\\t%AC_Het\\t%AC_Hom\\t%AC_Hemi\\n'
            cmd='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools +split-vep -d -f "{fields}"'.format(input_file, path=pop,fields=split_fields,file=vcf)
            cmd2='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools +fill-tags | /share/apps/bio/bin/bcftools query -f "{que}"'.format(input_file, path=pop, file=vcf, que=query)
            o1=sub.check_output(cmd, shell=True, encoding='utf8')
            o2=sub.check_output(cmd2, shell=True, encoding='utf8')
            with open("split_vep.txt", "a+") as out, open("fill_tags.txt", "a+") as out2:
                out.write(o1)
                out2.write(o2)
                out.close()
                out2.close()
        else:
            split_fields="%ID\\t%SYMBOL\\t%Feature\\t%Consequence\\t%AF\\t%CDS_position\\t%EXON\\t%INTRON\\t%Protein_position\\t%Amino_acids\\t%SIFT\\t%PolyPhen\\t%CADD_RAW[\\t%SAMPLE=%GT]\\n"
            awk_carr="""awk 'BEGIN{FS="\\t"}{for(i=1;i<14;i++) {printf "%s ",$i };for(j=14;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",$j};printf "\\n";}' """
            tr='tr " " "\t"'
            query='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\t%MAF\\t%AC_Het\\t%AC_Hom\\t%AC_Hemi\\n'
            cmd='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools +split-vep -d -f "{fields}" | {awk} | {tran}'.format(input_file, path=pop,fields=split_fields,file=vcf, awk=awk_carr, tran=tr)
            cmd2='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools +fill-tags | /share/apps/bio/bin/bcftools query -f "{que}"'.format(input_file, path=pop, file=vcf, que=query)
            o1=sub.check_output(cmd, shell=True, encoding='utf8')
            o2=sub.check_output(cmd2, shell=True, encoding='utf8')
            with open("split_vep.txt", "a+") as out, open("fill_tags.txt", "a+") as out2:
                out.write(o1)
                out2.write(o2)
                out.close()
                out2.close()
        return

    def sample_grepper(file):
        newli="\n"
        cmd="bcftools query -l {infile}".format(infile=file)
        out=sub.check_output(cmd, shell=True, encoding='utf8')
        sampl=out.split(newli)[:-1]
        return sampl

    def merge_and_filter(fill, s_vep):
        lis=sample_grepper(pop + "/" + "1.vcf.gz")
        fst_hea=["chr", "ps", "id", "ref", "alt", "ac", "an", "af", "maf", "ac_het", "ac_hom", "ac_hemi"]
        left=pd.read_table(fill, sep="\t", names=fst_hea,engine='c')
        snd_hea=["id", "sym", "feat", "cons", "gnomAF", "CDSps", "exon", "intr", "prot_ps", "AmAc", "sift", "polyphen", "cadd"]
        snd_hea.extend(lis)
        right=pd.read_table(s_vep, sep="\t", names=snd_hea,engine='c')
        merge=pd.merge(left, right, on='id', how="inner")
        dup=merge.drop_duplicates(subset=["chr", "ps", "id", "ref", "alt", "cons"])
    #dup=merge.drop_duplicates()
        effect_var=["splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "stop_lost", "start_lost", "stop_gained", "missense_variant", "NMD_transcript_variant", "protein_altering_variant", "splice_region_variant"]
        if dup["cons"].isin(effect_var).any():
            with open("Variant_search_" + pop_input + ".txt", "a+") as final:
                pd.set_option('max_colwidth', -1)
                out=dup[dup["cons"].isin(effect_var)]
                #drop=out[out["carr"].notna()]
                head=["CHROM", "POS", "ID", "REF", "ALT", "AC", "AN", "AF", "MAF", "AC_Het", "AC_Hom", "AC_Hemi", "Gene", "Isoform", "Consequence", "gnomAD_AF", "CDS_position", "Exon", "Intron", "Protein_position", "Amino_acids", "SIFT", "PolyPhen", "CADD_RAW"]
                head.extend(lis)
                data_ok=out.to_string(header=head, na_rep="NA", index=False)
                final.write(data_ok)
                final.close()
        else:
            print("\nThe research doesn't contain variants annotated as:")
            print("\n" + ', '.join(effect_var))
            dec=input(str("Do you want to print file anyway?\n\nDigit y/n: " ))
            if dec == "y":
                with open("Variant_search_" + pop_input + ".txt", "a+") as final:
                    pd.set_option('max_colwidth', -1)
                    drop=dup[dup["carr"].notna()]
                    head=["CHROM", "POS", "ID", "REF", "ALT", "AC", "AN", "AF", "MAF", "AC_Het", "AC_Hom", "AC_Hemi", "Gene", "Isoform", "Consequence", "gnomAD_AF", "CDS_position", "Exon", "Intron", "Protein_position", "Amino_acids", "SIFT", "PolyPhen", "CADD_RAW", "Carriers"]
                    lis=sample_grepper(pop + "/" + "1.vcf.gz")
                    sams=head + lis
                    data_ok=drop.to_string(header=sams, na_rep="NA", index=False)
                    final.write(data_ok)
                    final.close()
            else:
                os.remove(fill)
                os.remove(s_vep)
                exit(1)
        return

#def uniq_lines(infile):
    #lines_seen = set() # holds lines already seen
    #outfile = open("Variant_search_" + pop_input + ".txt", "w")
    #for line in open(infile, "r"):
        #if line not in lines_seen: # not a duplicate
            #outfile.write(line)
            #lines_seen.add(line)
    #outfile.close()
    #return

    def imputed_search(vcf):
        if pop_input != "CAR" and pop_input != "VBI" and pop_input != "I_FVG":
            try:
                query='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/IMP\\t%INFO/AF\\t%INFO/INFO[\\t%SAMPLE=%DS]\\n'
                awk_carr="""awk 'BEGIN{FS="\\t"}{for(i=1;i<9;i++) {printf "%s ",$i };for(j=9;j<=NF;j++){split($j,a,"=");if ((a[2] >= 0.9)) printf "%s,",$j};printf "\\n";}' """
                tr='tr " " "\t"'
                cmd='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools query -f "{que}" | {awk} | {tran}'.format(input_file, path=pop, awk=awk_carr, tran=tr, file=vcf, que=query)
                o1=sub.check_output(cmd, shell=True, encoding='utf8')
            except:
                query='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/IMP\\t%INFO/AF\\t%INFO/R2[\\t%SAMPLE=%DS]\\n'
                awk_carr="""awk 'BEGIN{FS="\\t"}{for(i=1;i<9;i++) {printf "%s ",$i };for(j=9;j<=NF;j++){split($j,a,"=");if ((a[2] >= 0.9)) printf "%s,",$j};printf "\\n";}' """
                tr='tr " " "\t"'
                cmd='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools query -f "{que}" | {awk} | {tran}'.format(input_file, path=pop, awk=awk_carr, tran=tr, file=vcf, que=query)
                o1=sub.check_output(cmd, shell=True, encoding='utf8')
            else:
                query='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/IMP2_exp_freq_a1\\t%INFO/R2[\\t%SAMPLE=%DS]\\n'
                awk_carr="""awk 'BEGIN{FS="\\t"}{for(i=1;i<8;i++) {printf "%s ",$i };for(j=8;j<=NF;j++){split($j,a,"=");if ((a[2] >= 0.9)) printf "%s,",$j};printf "\\n";}' """
                tr='tr " " "\t"'
                cmd='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools query -f "{que}" | {awk} | {tran}'.format(input_file, path=pop, awk=awk_carr, tran=tr, file=vcf, que=query)
                o1=sub.check_output(cmd, shell=True, encoding='utf8')
                with open("imputed_search_" + pop_input + ".txt", "a") as f_out:
                    f_out.write(o1)
                    f_out.close()
                    return


#if build_input != 37 and build_input != 38 or data_type != "WGS" and data_type != "IMP":
    #print("***ERROR***\n\nPlease select a valid build reference and type of data.\nBuild reference: 37 or 38\nData type: WGS or IMP")
    #exit(1)


    if build_input == 37:
        data_type=str(input("Please insert the type of data you want to query.\n\nChoose between WGS or IMP (imputed data): "))
        if data_type == "WGS":
            pop_input=str(input("Please insert the desired population.\nChoose between:\nFVG, CAR, VBI: "))
            pop=pop_path37[pop_input]
            crom=list(regions.cro)
            for chr in crom:
                vcf=chr + ".vcf.gz"
                bcf_search(vcf)
                merge_and_filter("fill_tags.txt", "split_vep.txt")
        #uniq_lines("var_screen_" + pop_input + ".txt")
            os.remove("fill_tags.txt")
            os.remove("split_vep.txt")
        else:
            pop_input=str(input("\nPlease insert the desired population from imputed data.\n\nPlease choose from MOL, ITT, CAR, VBI, I_FVG, DIABETES, SLO, NUTRIACT, ERB, VBI, BAU: "))
            pop=imputed_data[pop_input]
            head='CHR\tPOS\tID\tREF\tALT\tImputed_Marker\tEstimated_AF\tInfo_score\tGenotype_Dosage'
            x=0
            while x < 22:
                if pop_input == "MOL" or pop_input == "ITT":
                    chr=str(CHR[x])
                    vcf="/" + chr + "/" + chr + ".vcf.gz"
                    imputed_search(vcf)
                    x=x+1
                elif pop_input == "CAR" or pop_input =="VBI" or pop_input == "I_FVG":
                    chr=str(CHR[x])
                    vcf="/chr" + chr + ".vcf.gz"
                    head='CHR\tPOS\tID\tREF\tALT\tEstimated_AF\tInfo_score\tGenotype_Dosage'
            #print(pop+vcf)
                    imputed_search(vcf)
                    x=x+1
                else:
                    chr=str(CHR[x])
                    vcf="/" + chr + ".vcf.gz"
                    imputed_search(vcf)
                    x=x+1
                    with open("imputed_search_" + pop_input + ".txt", "r") as out:
                        save=out.read()
                        out.close()
                        with open("imputed_search_" + pop_input + ".txt", "w") as new_out:
                            new_out.write(head+"\n")
                            new_out.write(save)
                            new_out.close()
    else:
        pop_input=str(input("\nOnly WGS data available. Please insert the desired population.\nChoose between:\nSR, HCFVG, COV, OTO, WCS: "))
        pop=pop_path38[pop_input]
        crom=list(regions.cro)
        for chr in crom:
            final_chr=chr[3:]
            vcf=final_chr + ".vcf.gz"
            bcf_search(vcf)
            merge_and_filter("fill_tags.txt", "split_vep.txt")
        os.remove("fill_tags.txt")
        os.remove("split_vep.txt")

#if sys.stderr.getvalue():
   #with open (pop_input + ".e") as output:
        #output.write(sys.stderr.getvalue())
        #output.close()

#JOB mode

elif mode == "JOB":
    cbuild=sys.argv[3]
    def bcf_search(vcf):
        if  compop != "SR":
            split_fields="%CHROM\\t%POS\\t%ID\\t%QUAL\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\t%INFO/AF\\t%INFO/MAF\\t%INFO/AC_Het\\t%INFO/AC_Hom\\t%INFO/AC_Hemi\\t%SYMBOL\\t%Feature\\t%Consequence\\t%gnomAD_AF\\t%CDS_position\\t%EXON\\t%INTRON\\t%Protein_position\\t%Amino_acids\\t%SIFT\\t%PolyPhen\\t%CADD_RAW[\\t%GT]\\n"
            #awk_carr="""awk 'BEGIN{FS="\\t"}{for(i=1;i<14;i++) {printf "%s ",$i };for(j=14;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",$j};printf "\\n";}' """
            #tr='tr " " "\t"'
            query='%ID\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\t%MAF\\t%AC_Het\\t%AC_Hom\\t%AC_Hemi\\n'
            cmd='/share/apps/bio/bin/bcftools view -i "AC_Hom != 0" -r {chr}:{start}-{end} {path}/{file} | /share/apps/bio/bin/bcftools +split-vep -d -f "{fields}"'.format(chr=chr,start=start,end=end, path=pop,fields=split_fields,file=vcf)
            #cmd2='/share/apps/bio/bin/bcftools view -r {chr}:{start}-{end} {path}/{file} | /share/apps/bio/bin/bcftools +fill-tags | /share/apps/bio/bin/bcftools query -f "{que}"'.format(chr=chr,start=start,end=end, path=pop, file=vcf, que=query)
            o1=sub.check_output(cmd,shell=True, encoding='utf8')
            #o2=sub.check_output(cmd2,shell=True, encoding='utf8')
            with open("/netapp05/analisi_nardone/COV19/split_vep.txt", "a+") as out:
                out.write(o1)
                #out2.write(o2)
                out.close()
                #out2.close()
        else:
            split_fields="%CHROM\\t%ID\\t%SYMBOL\\t%Feature\\t%Consequence\\t%AF\\t%CDS_position\\t%EXON\\t%INTRON\\t%Protein_position\\t%Amino_acids\\t%SIFT\\t%PolyPhen\\t%CADD_RAW[\\t%GT]\\n"
            #awk_carr='BEGIN{FS="\\t"}{for(i=1;i<14;i++) {printf "%s ",$i };for(j=14;j<=NF;j++){split($j,a,"=");if ((a[2]!="0|0" && a[2]!=".|.") && (a[2]!="0/0" && a[2]!="./.")) printf "%s,",$j};printf "\\n";}' """
            #tr='tr " " "\t"'
            query='%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AC\\t%AN\\t%AF\\t%MAF\\t%AC_Het\\t%AC_Hom\\t%AC_Hemi\\n'
            cmd='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools +split-vep -d -f "{fields}"'.format(input_file, path=pop,fields=split_fields,file=vcf)
            cmd2='/share/apps/bio/bin/bcftools view -R {} {path}/{file} | /share/apps/bio/bin/bcftools +fill-tags | /share/apps/bio/bin/bcftools query -f "{que}"'.format(input_file, path=pop, file=vcf, que=query)
            o1=sub.check_output(cmd, shell=True, encoding='utf8')
            o2=sub.check_output(cmd2, shell=True, encoding='utf8')
            with open("split_vep.txt", "a+") as out, open("fill_tags.txt", "a+") as out2:
                out.write(o1)
                out2.write(o2)
                out.close()
                out2.close()
        return

    def sample_grepper(file):
        newli="\n"
        cmd="bcftools query -l {infile}".format(infile=file)
        out=sub.check_output(cmd, shell=True, encoding='utf8')
        sampl=out.split(newli)[:-1]
        return sampl

    def merge_and_filter(fill, s_vep):
        lis=sample_grepper(pop + "/" + "1.vcf.gz")
        fst_hea=["chr","ps","id","qual","ref", "alt", "ac", "an", "af", "maf", "ac_het", "ac_hom", "ac_hemi","sym", "feat", "cons", "gnomAF", "CDSps", "exon", "intr", "prot_ps", "AmAc", "sift", "polyphen", "cadd"]
        #left=pd.read_csv(fill, sep="\t", names=fst_hea,engine='c',low_memory=False,dtype='object')
        #snd_hea=["sym", "feat", "cons", "gnomAF", "CDSps", "exon", "intr", "prot_ps", "AmAc", "sift", "polyphen", "cadd"]
        fst_hea.extend(lis)
        right=pd.read_csv(s_vep, sep="\t", names=fst_hea,engine='c',dtype='object',low_memory=False)
        right1=right.compute()
        effect_var=["splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "stop_lost", "start_lost", "stop_gained","NMD_transcript_variant", "protein_altering_variant", "splice_region_variant"]
        if right1["cons"].isin(effect_var).any():
            lof=right1.query("cons in @effect_var")
            miss=right1.query("cons=='missense_variant'")
            #lof_merge=.drop_duplicates(subset=['id','ps']).compute()
            #miss_merge.drop_duplicates(subset=['id','ps']).compute()
            head=["CHROM", "POS", "ID","PHRED","REF", "ALT", "AC", "AN", "AF", "MAF", "AC_Het", "AC_Hom", "AC_Hemi", "Gene", "Isoform", "Consequence", "gnomAD_AF", "CDS_position", "Exon", "Intron", "Protein_position", "Amino_acids", "SIFT", "PolyPhen", "CADD_RAW"]
            head.extend(lis)
            with open("Variant_search_LoF_"+compop+pers_name+".txt",'a+') as lf, open("Variant_search_missense_"+compop+pers_name+".txt",'a+') as ms:
                out_lof=lof.to_csv(header=head,na_rep="NA", index=False, sep="\t")
                out_miss=miss.to_csv(header=head,na_rep="NA", index=False, sep="\t")
                lf.write(out_lof)
                ms.write(out_miss)
                lf.close()
                ms.close()
        #dup=merge.drop_duplicates(subset=["chr", "ps", "id", "ref", "alt", "cons"])
    #dup=merge.drop_duplicates()
        #effect_var=["splice_acceptor_variant", "splice_donor_variant", "frameshift_variant", "stop_lost", "start_lost", "stop_gained", "missense_variant", "NMD_transcript_variant", "protein_altering_variant", "splice_region_variant"]
        #dup1=merge.compute()
        #if dup1["cons"].isin(effect_var).any():
            #with open("Variant_search_" + compop + pers_name + ".txt", "a+") as final:
                #pd.set_option('max_colwidth', -1)
                #out=dup1[dup1["cons"].isin(effect_var)]
                #drop=out[out["carr"].notna()]
                #head=["CHROM", "POS", "ID", "REF", "ALT", "AC", "AN", "AF", "MAF", "AC_Het", "AC_Hom", "AC_Hemi", "Gene", "Isoform", "Consequence", "gnomAD_AF", "CDS_position", "Exon", "Intron", "Protein_position", "Amino_acids", "SIFT", "PolyPhen", "CADD_RAW"]
                #head.extend(lis)
                #data_ok=out.to_string(header=head, na_rep="NA", index=False)
                #final.write(data_ok)
                #final.close()
        else:
            print("\nThe research doesn't contain variants annotated as:")
            print("\n" + ', '.join(effect_var)+"or missense variant")
            print("\n"+"Printing relults anyway")
            with open("Variant_search_" + compop + pers_name + ".txt", "a+") as final:
                pd.set_option('max_colwidth', -1)
                #drop=dup[dup["carr"].notna()]
                head=["CHROM", "POS", "ID", "REF", "ALT", "AC", "AN", "AF", "MAF", "AC_Het", "AC_Hom", "AC_Hemi", "Gene", "Isoform", "Consequence", "gnomAD_AF", "CDS_position", "Exon", "Intron", "Protein_position", "Amino_acids", "SIFT", "PolyPhen", "CADD_RAW", "Carriers"]
                head.extend(lis)
                data_ok=out.to_string(header=head, na_rep="NA", index=False)
                final.write(data_ok)
                final.close()
        return
    if cbuild == "37":
        cdata=sys.argv[4]
        compop=sys.argv[5]
        pers_name=sys.argv[6]
        #data_type=str(input("Please insert the type of data you want to query.\n\nChoose between WGS or IMP (imputed data): "))
        if cdata == "WGS":
            #pop_input=str(input("Please insert the desired population.\nChoose between:\nFVG, CAR, VBI: "))
            pop=pop_path37[compop]
            crom=list(regions.cro)
            for chr in crom:
                vcf=chr + ".vcf.gz"
                bcf_search(vcf)
            merge_and_filter("fill_tags.txt", "split_vep.txt")
        #uniq_lines("var_screen_" + pop_input + ".txt")
            os.remove("fill_tags.txt")
            os.remove("split_vep.txt")
        else:
            #pop_input=str(input("\nPlease insert the desired population from imputed data.\n\nPlease choose from MOL, ITT, CAR, VBI, I_FVG, DIABETES, SLO, NUTRIACT, ERB, VBI, BAU: "))
            pop=imputed_data[compop]
            head='CHR\tPOS\tID\tREF\tALT\tImputed_Marker\tEstimated_AF\tInfo_score\tGenotype_Dosage'
            x=0
            while x < 22:
                if compop == "MOL" or compop == "ITT":
                    chr=str(CHR[x])
                    vcf="/" + chr + "/" + chr + ".vcf.gz"
                    imputed_search(vcf)
                    x=x+1
                elif compop == "CAR" or compop =="VBI" or compop == "I_FVG":
                    chr=str(CHR[x])
                    vcf="/chr" + chr + ".vcf.gz"
                    head='CHR\tPOS\tID\tREF\tALT\tEstimated_AF\tInfo_score\tGenotype_Dosage'
            #print(pop+vcf)
                    imputed_search(vcf)
                    x=x+1
                else:
                    chr=str(CHR[x])
                    vcf="/" + chr + ".vcf.gz"
                    imputed_search(vcf)
                    x=x+1
                    with open("imputed_search_" + compop + ".txt", "r") as out:
                        save=out.read()
                        out.close()
                        with open("imputed_search_" + compop + ".txt", "w") as new_out:
                            new_out.write(head+"\n")
                            new_out.write(save)
                            new_out.close()
    elif cbuild == "38":
        #pop_input=str(input("\nOnly WGS data available. Please insert the desired population.\nChoose between:\nSR, HCFVG, COV, OTO, WCS: "))
        compop=sys.argv[4]
        pers_name=sys.argv[5]
        pop=pop_path38[compop]
        #crom=list(regions.cro)
        #for chr in crom:
            #final_chr=chr[3:]
            #vcf=final_chr + ".vcf.gz"
        if compop != "HCFVG2":
            for line in open(input_file,'r'):
                chr=line.split("\t")[0]
                start=line.split("\t")[1]
                end=line.split("\t")[2].replace("\n","")
                final_chr=chr[3:]
                vcf=final_chr + ".vcf.gz"
                bcf_search(vcf)
        else:
            for line in open(input_file,'r'):
                chr=line.split("\t")[0]
                start=line.split("\t")[1]
                end=line.split("\t")[2].replace("\n","")
                vcf=chr + ".vcf.gz"
                bcf_search(vcf)
        merge_and_filter("fill_tags.txt", "split_vep.txt")
        #os.remove("fill_tags.txt")
        os.remove("split_vep.txt")

#if sys.stderr.getvalue():
   #with open (pop_input + ".e") as output:
        #output.write(sys.stderr.getvalue())
        #output.close()
