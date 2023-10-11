import csv
import sys
import random
import os
from collections import Counter
from numpy.random import randn

target_gene_id=sys.argv[1]
num_rare=int(sys.argv[2])
num_major=5


frank_region = 500000

with open("/PATH/TO/gencode.v26.GRCh38.csv") as gencode_csv: ## Path to genocode file
    readgencodeCSV = csv.reader(gencode_csv, delimiter=',')
    for row in readgencodeCSV:
        curr_gene_id = row[3].strip()
        if curr_gene_id == target_gene_id:
            target_chr=row[0].strip()
            loc_start=int(row[1].strip())
            loc_end=int(row[2].strip())
            break
gencode_csv.close()

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
        if loc_start <= int(curr_tmp_loc) <=loc_end:
            rare_var_list.append(line.split("\t")[0].strip()+"_"+line.split("\t")[3].strip())
rare_var_file.close()

overlap_file=open("/PATH/TO/HG_GTEx_v8_chr"+target_chr+".bim") ## The bim file for overlapping variants between GTEx and 1000HG
snp_b38_list=[]
chr_list=[]
position_list=[]
chr_pos_list=[]
nonrare_var_list=[]
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
            if curr_chr_pos not in rare_var_list:
                nonrare_var_list.append(curr_chr_pos)
            chr_pos_list.append(curr_chr_pos)
            snp_b38_list.append(line.split("\t")[1].strip())
            chr_list.append(curr_tmp_chr)
            position_list.append(curr_tmp_loc)
            chrpos_a0_dict[curr_chr_pos]=line.split("\t")[5].strip()
            chrpos_a1_dict[curr_chr_pos]=line.split("\t")[4].strip()
            chrpos_rs_dict[curr_chr_pos]=line.split("\t")[1].strip()
            if curr_chr_pos in rare_var_list:
                chrpos_maf_dict[curr_chr_pos]=str(0.005)
            else:
                chrpos_maf_dict[curr_chr_pos]=str(0.05)


random_index_rare=random.sample(range(len(rare_var_list)),num_rare)
num_non_rare=num_major
print(num_non_rare)
print(num_rare)
rare_chrpos_list=[]
for index in range(len(rare_var_list)):
    if index in random_index_rare:
        rare_chrpos_list.append(rare_var_list[index])
print("rare_list")
print(len(rare_chrpos_list))

random_index_nonrare=random.sample(range(len(nonrare_var_list)),num_non_rare+1)
random_SEL_list=random_index_nonrare[1:]
random_COT=random_index_nonrare[0]
nonrare_chrpos_list=[]
for index in range(len(nonrare_var_list)):
    if index in random_SEL_list:
        nonrare_chrpos_list.append(nonrare_var_list[index])
    if index ==random_COT:
        COT_chrpos=nonrare_var_list[index]
print("nonrare_list")
print(len(nonrare_chrpos_list))

# Generatio COEF from normal distribution
#print(num_causal)
#random_COF_SEL_list=randn(num_causal)
#random_COF_COT=randn(1)

out_dir="/PATH/TO/OUT_DIR/"+target_gene_id+"/" ##Path to the output directory from the first step
output=open(out_dir+"INFO_FUNC_"+str(num_rare)+"_"+target_gene_id+"_ADD.txt",'w')
output.write("#CHR\tPOS\tRS_ID\tREF\tALT\tORIGIN\tCOEFF\tMAF\tSNP_b38")
tmp_chrpos_list=[]
for chrpos in chr_pos_list:
    if chrpos not in tmp_chrpos_list:
        output.write("\n"+chrpos.split("_")[0].strip())
        output.write("\t"+chrpos.split("_")[1].strip())
        output.write("\t"+chrpos_rs_dict[chrpos])
        output.write("\t"+chrpos_a0_dict[chrpos])
        output.write("\t"+chrpos_a1_dict[chrpos])
        if chrpos in nonrare_chrpos_list:
            output.write("\tSEL\t"+str(randn(1)[0])+"\t"+chrpos_maf_dict[chrpos])
        elif chrpos in rare_chrpos_list:
            output.write("\tSEL\t"+str(randn(1)[0])+"\t"+chrpos_maf_dict[chrpos])
        elif chrpos == COT_chrpos:
            output.write("\tCOT\t"+str(randn(1)[0])+"\t"+chrpos_maf_dict[chrpos])
        else:
            output.write("\tNULL\t0.0\t"+chrpos_maf_dict[chrpos])
        tmp_chrpos_list.append(chrpos)
output.close()

