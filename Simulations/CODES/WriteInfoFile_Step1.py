import csv
import sys
import random
import os
from collections import Counter
from numpy.random import randn

target_gene_id=sys.argv[1]

frank_region = 500000

with open("/PATH/TO/gencode.v26.GRCh38.csv") as gencode_csv:
    readgencodeCSV = csv.reader(gencode_csv, delimiter=',')
    for row in readgencodeCSV:
        curr_gene_id = row[3].strip()
        if curr_gene_id == target_gene_id:
            target_chr=row[0].strip()
            loc_start=int(row[1].strip())
            loc_end=int(row[2].strip())
            break
gencode_csv.close()
print(target_gene_id)
print(target_chr)

if loc_start < frank_region:
    loc_start = 1;
else:
    loc_start = loc_start - frank_region
loc_end = loc_end + frank_region

rare_var_file=open("/PATH/TO/HG_GTEx_v8_chr"+target_chr+"_rare.bim") ## The bim file for overlapping rare variants between GTEx and 1000HG
rare_var_list=[]
for line in rare_var_file:
    curr_tmp_chr = line.split("\t")[0].strip()
    if curr_tmp_chr == target_chr:
        curr_tmp_loc = line.split("\t")[3].strip()
        curr_chr_pos=curr_tmp_chr+"_"+curr_tmp_loc
        if loc_start <= int(curr_tmp_loc) <=loc_end:
            rare_var_list.append(curr_chr_pos)
rare_var_file.close()
print("rare_var_list")
print(len(rare_var_list))

overlap_file=open("/PATH/TO/HG_GTEx_v8_chr"+target_chr+".bim") ## The bim file for overlapping variants between GTEx and 1000HG
snp_b38_list=[]
chr_list=[]
position_list=[]
chr_pos_list=[]
chrpos_a0_dict={}
chrpos_a1_dict={}
chrpos_rs_dict={}
chrpos_maf_dict={}
for line in overlap_file:
    curr_tmp_chr = line.split("\t")[0].strip()
    if curr_tmp_chr == target_chr:
        curr_tmp_loc = line.split("\t")[3].strip()
        if loc_start <= int(curr_tmp_loc) <=loc_end:
            curr_chr_pos=curr_tmp_chr+"_"+curr_tmp_loc
            chr_pos_list.append(curr_chr_pos)
            snp_b38_list.append(line.split("\t")[1].strip())
            chr_list.append(curr_tmp_chr)
            position_list.append(curr_tmp_loc)
            chrpos_a0_dict[curr_chr_pos]=line.split("\t")[5].strip()
            chrpos_a1_dict[curr_chr_pos]=line.split("\t")[4].strip()
            chrpos_rs_dict[curr_chr_pos]=line.split("\t")[1].strip()
            if line.split("\t")[1].strip() in rare_var_list:
                chrpos_maf_dict[curr_chr_pos]=str(0.005)
            else:
                chrpos_maf_dict[curr_chr_pos]=str(0.05)

print("overlap")
print(len(chr_pos_list))

out_dir="/PATH/TO/OUT_DIR/"+target_gene_id+"/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

output_csv=open(out_dir+target_gene_id+"_gtex_genotype.csv",'w')
output_csv_maf=open(out_dir+target_gene_id+"_gtex_genotype_maf0.05.csv",'w')
gtex_orig_file="/PATH/TO/GTEX_Genotype/GTEx_v8_866Indiv.clean_chr"+target_chr+".num.csv" ##path to GTEx genotype file 
snp_info_dict={}
pheno_id_list=[]
with open(gtex_orig_file) as gtex_csv:
    readGTEX = csv.reader(gtex_csv, delimiter=',')
    header=next(readGTEX)
    output_csv.write("CHR,POS")
    output_csv_maf.write("CHR,POS")
    for idv in header[2:]:
        output_csv.write(","+idv.strip())
        output_csv_maf.write(","+idv.strip())
        pheno_id_list.append(idv.strip())
    for row in readGTEX:
        curr_b38=row[0].strip()+"_"+row[1].strip()
        if curr_b38 in chr_pos_list:
            curr_snp_info=row[0].strip()+","+row[1].strip()
            for curr_snp in row[2:]:
                curr_snp_info=curr_snp_info+","+curr_snp.strip()
            snp_info_dict[curr_b38]=curr_snp_info

tmp_chrpos_list=[]
for curr_snp in chr_pos_list:
    if curr_snp not in tmp_chrpos_list:
        output_csv.write("\n"+snp_info_dict[curr_snp])
        if curr_snp not in rare_var_list:
            output_csv_maf.write("\n"+snp_info_dict[curr_snp])
        tmp_chrpos_list.append(curr_snp)
output_csv.close()
gtex_csv.close()

hg_csv_file="/PATH/TO/1000HG_Genotype/ALL.chr"+target_chr+".GRCh38.phased.clean.num.csv" ##path to 1000HG genotype file
output_hgtped=open(out_dir+target_gene_id+"_hg1000_dosage.txt",'w')
output_hgcsv=open(out_dir+target_gene_id+"_hg1000_genotype.csv",'w')
output_hgcsv_maf=open(out_dir+target_gene_id+"_hg1000_genotype_maf0.05.csv",'w')
output_hgcsv_mafnew=open(out_dir+target_gene_id+"_hg1000_genotype_maf0.05_new.csv",'w')
hg_pheno_id=[]
hgsnp_info_dict={}
hgtped_info_dict={}
with open(hg_csv_file) as hg_csv:
    readHG = csv.reader(hg_csv, delimiter=',')
    header=next(readHG)
    output_hgcsv.write("CHR,POS")
    output_hgcsv_maf.write("CHR,POS")
    output_hgcsv_mafnew.write("CHR,POS")
    for idv in header[2:]:
        output_hgcsv.write(","+idv.strip())
        output_hgcsv_maf.write(","+idv.strip())
        hg_pheno_id.append(idv.strip())
    for row in readHG:
        curr_b38=row[0].strip()+"_"+row[1].strip()
        if curr_b38 in chr_pos_list:
            curr_chr_pos=curr_b38
            curr_chr=row[0].strip()
            curr_pos=row[1].strip()
            curr_a0=chrpos_a0_dict[curr_chr_pos]
            curr_a1=chrpos_a1_dict[curr_chr_pos]
            curr_rs=chrpos_rs_dict[curr_chr_pos]
            curr_maf=chrpos_maf_dict[curr_chr_pos]
            hg_snp_tped=curr_chr+" "+curr_rs+" "+curr_pos+" "+curr_a0+" "+curr_a1+" "+curr_maf
            curr_snp_info=row[0].strip()+","+row[1].strip()
            for curr_char in row[2:]:
                curr_snp_info=curr_snp_info+","+curr_char.strip()
                hg_snp_tped=hg_snp_tped+" "+curr_char.strip()
            hgsnp_info_dict[curr_b38]=curr_snp_info
            hgtped_info_dict[curr_b38]=hg_snp_tped

tmp_chrpos_list=[]
for curr_snp in chr_pos_list:
    if curr_snp not in tmp_chrpos_list:
        output_hgcsv.write("\n"+hgsnp_info_dict[curr_snp])
        output_hgtped.write(hgtped_info_dict[curr_snp]+"\n")
        if curr_snp not in rare_var_list:
            output_hgcsv_maf.write("\n"+hgsnp_info_dict[curr_snp])
            output_hgcsv_mafnew.write("\nchr"+hgsnp_info_dict[curr_snp])
        tmp_chrpos_list.append(curr_snp)
output_hgcsv.close()
hg_csv.close()
