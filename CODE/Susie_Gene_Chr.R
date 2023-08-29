library(susieR)
library(methods)
library(glmnet)
library(stringr)
library(data.table)
library(dplyr)

"%&%" <- function(a,b) paste(a,b, sep = "")

args = commandArgs(trailingOnly=TRUE)

chrom<-args[1]
prefix<-"GTEx_Model"

#### I. data file preparation #### 
snp_annot_file <- '/PATH/TO/GTEx_v8_866Indiv_wb.clean_chr'%&% chrom %&%'.bim'
gene_annot_file <- "/PATH/TO/gencode.v26.GRCh38.genes.gtf"
genotype_file <- "/PATH/TO/GTEx_v8_866Indiv_wb.clean_chr"%&% chrom %&% ".num.csv" 
expression_file <- "/PATH/TO//Whole_Blood.v8.residuals.csv" ## gene epxression bed file
gwas_file <- "/path/to/gwas_file/"

#### II. load data  ######
# snp_annot_file from GTEx-those variants that rank as top 50k
# Make different snp_anno_file for 50k,100k...all_snp, contain all the overlapping variants between GTEx and GWAS that column with
snp_annot <- read.table(snp_annot_file, header = T, stringsAsFactors = F,sep="\t")
colnames(snp_annot)=c("chr","varID","0","pos","alt","ref")
snp_annot$chrpos<-snp_annot$chr %&%"_"%&% snp_annot$pos
gwassnp <- fread(file=gwas_file,header=T,sep=",")
gwassnp <- as.data.frame(gwassnp)
gwas_snp<- gwassnp$CHR%&%"_"%&%gwassnp$LOC
snp_annot<-snp_annot[snp_annot$chrpos %in% gwas_snp,]
print(dim(snp_annot))

get_gene_annotation <- function(gene_annot_file_name, chrom)
{
  gene_df <- read.table(gene_annot_file,header=F,stringsAsFactors =F,sep="\t",fill=T)
  gene_df1 <- filter(gene_df,V3 %in% "gene")
  geneid <- str_extract(gene_df1[,9], "ENSG\\d+.\\d+")
  genename <- gsub("gene_name (\\S+);","\\1",str_extract(gene_df1[,9], "gene_name (\\S+);"), perl=T)
  gene_used <- as.data.frame(cbind(geneid,genename,gene_df1[,c(1,4,5,3)]))
  colnames(gene_used) <- c("geneid","genename","chr","start","end","anno")
  gtf_used <- filter(gene_used,gene_used[,3] %in% ('chr' %&% chrom))
  gtf_used
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  #expr_df <- as.data.frame(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1))
  expr_df <- as.data.frame(read.table(file=gene_expression_file_name,header=FALSE,sep=",",check.names = FALSE))
  col_name=expr_df$V4
  expr_df<- expr_df[,c(-1,-2,-3,-4)]
  expr_df_T = as.data.frame(t(expr_df))
  colnames(expr_df_T)=col_name
  rownames(expr_df_T)=colnames(expr_df)
  expr_df_T <- expr_df_T %>% select(one_of(intersect(gene_annot$geneid, colnames(expr_df_T))))
  expr_df_T
}

get_filtered_genotype <- function(genotype_file_name) {
  gt_df <- fread(file=genotype_file_name,header=T,sep=",")
  gt_df <- as.data.frame(gt_df)
  snp_name<-gt_df$CHR %&% "_"%&%gt_df$LOC
  gt_df_T<- as.data.frame(t(as.matrix(gt_df[, c(-1,-2)])))
  colnames(gt_df_T) <- snp_name 
  gt_df_T
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$geneid == gene),]
  c(row$start, row$end)
}


get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window)  & (pos <= (coords[2] + cis_window))))
  #print(colnames(gt_df))
  cis_gt <- gt_df[, colnames(gt_df) %in% snp_info$chrpos]
  cis_gt
}

### main function #####
cis_window=500000
# Read in data----
gene_annot <- get_gene_annotation(gene_annot_file, chrom)
print(dim(gene_annot))
expr_df <- get_gene_expression(expression_file, gene_annot)
print(dim(expr_df))
samples <- rownames(expr_df)
n_samples <- length(samples)
#print(n_samples)
genes <- colnames(expr_df)
print(genes[1:10])
n_genes <- length(expr_df)
print(n_genes)
gt_df <- get_filtered_genotype(genotype_file)
print(dim(gt_df))

# Prepare output data----
weights_file <- './summary/' %&% prefix %&% '_chr' %&% chrom %&% '_weights.txt'
weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
write(weights_col, file = weights_file, ncol = 6, sep = '\t')


for (i in 1:n_genes) {
  cat(i, "/", n_genes, "\n")
  gene <- genes[i]
  print(gene)
  gene_name <- gene_annot$genename[gene_annot$geneid == gene]
  #print(gene_name)
  coords <- get_gene_coords(gene_annot, gene)
  cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)
  if(length(dim(cis_gt))==2){
  	if (ncol(cis_gt) >= 2) {
		for(k in 1:nrow(cis_gt)) {       # for-loop over columns
        		cis_gt[k ,] <- as.double(cis_gt[k , ])
  		}
  		cis_gt=data.matrix(cis_gt)
  		print(class(cis_gt))
  		print(dim(cis_gt))
  		print(typeof(cis_gt))
  		snp_name_list<-colnames(cis_gt)
    	adj_expression <- expr_df[,i]
		adj_expression=as.vector(adj_expression)
		X=cis_gt
        y=adj_expression
        fitted = susieR::susie(X, y, L=min(10,ncol(X)))
		bhat = coef(fitted)[-1]
		snp_pip_list=fitted$pip
		print(fitted$sets$cs)
		print(fitted$sets$cs$L1)
		index=1
                sel_weight_list=c()
                sel_snp_list=c()
		sel_pip_list=c()
		if (length(fitted$sets$cs) != 0){
                        set_len=length(fitted$sets$cs)
                        for (ss in 1:set_len){
                                list_we=fitted$sets$cs[[ss]]
                                for (v in 1:length(snp_name_list)) {
                                        if (v %in% list_we){
						if (!snp_name_list[v] %in% sel_snp_list){
                                                	sel_snp_list[index]=snp_name_list[v]
                                                	sel_weight_list[index]=bhat[v]
                                                	sel_pip_list[index]=snp_pip_list[v]
                                                	index=index+1
						}
                                        }
                                }
                        }
                }
      		if (length(sel_weight_list) > 0) {
        		weighted_snps_info <- snp_annot %>% filter(chrpos %in% sel_snp_list) %>% select(varID,chrpos,ref,alt)
        		print("weighted_snp_info")
			print(weighted_snps_info)
			if (nrow(weighted_snps_info) != 0){
        			weighted_snps_info$gene <- gene
        			weighted_snps_info <- weighted_snps_info %>% merge(data.frame(weights =sel_weight_list, chrpos=sel_snp_list), by = 'chrpos') %>%
          			select(gene, varID, chrpos, ref, alt, weights)
        			write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
			}
      		}
    	}
  }
}

