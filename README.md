## Users’ Manual of rvTWAS
### Overview:
Towards the identification of genetic basis of complex traits, transcriptome-wide association study (TWAS) is successful in integrating transcriptome data. However, TWAS is only applicable for common variants, excluding rare variants in exome or whole genome sequences. This is partly because of the inherent limitation of TWAS protocols that rely on predicting gene expressions. Briefly, a typical TWAS protocol has two steps: it trains an expression prediction model in a reference dataset containing gene expressions and genotype, and then applies this prediction model to a genotype-phenotype dataset to “impute” the unobserved expression (that is called GReX) to be associated to the phenotype. In this procedure, rare variants are not used due to its low power in predicting expressions. Our previous research has revealed the insight into TWAS: the two steps are essentially genetic feature selection and aggregations that do not have to involve predictions. Based on this insight disentangling TWAS, rare variants’ inability of predicting expression traits is no longer an obstacle. Herein, we developed “rare variant TWAS”, or rvTWAS, that first uses a Bayesian model to conduct expression-directed feature selection and then use a kernel machine to carry out feature aggregation, forming a model leveraging expressions for association mapping including rare variants. We demonstrated the performance of rvTWAS by thorough simulations and real data analysis in three psychiatric disorders, namely schizophrenia, bipolar disorder, and autism spectrum disorder. rvTWAS will open a door for sequence-based association mappings integrating gene expressions.

![My Image](EXAMPLE/Fig1A_B.PNG)

### Installation
rvTWAS is a batteries-included JAR executable. All needed external jar packages are included in the downloadable, rvTWAS.jar. Then users can use the commands:\
`git clone https://github.com/QingrunZhangLab/rvTWAS.git`

As we used an R package "susieR" and "SKAT", the users have to install "R", "susieR"(https://cran.r-project.org/web/packages/susieR/susieR.pdf),and "SKAT"(https://cran.r-project.org/web/packages/SKAT/index.html). The versions of "R", "SuSiE", and "SKAT" packages that we have used on our platform are: version 3.6.1 for "R", version 0.12.35 for "susieR", and version 2.0.0 for "SKAT". Users are also expected to have java (version: 1.8) and Plink (version: 1.9) installed on their platform.

Usage:
rvTWAS, which is composed by two steps: First, rvTWAS uses SuSiE[1] to carry out variants selections to form a prioritized set of genetic variants (including rare variants) weighted by their relevance to gene expressions (Figure 1A). Second, supported by our previous  successful attempt using kernel methods to carry out common variants-based TWAS, as well as the practice of using kernel models in both common and rare variants GWAS, we use SKAT[2][3] method to aggregate weighted variants to form a score test for the association (Figure 1B). 

### 1. Feature selection using SuSiE:
#### 1.1 Prepare input data:
**1.1.1	Gene expression file:** \
Take the whole blood as an example. The fully processed, filtered and normalized gene expression matrices in bed format ("Whole_Blood.v8.normalized_expression.bed") for whole blood was downloaded from GTEx portal (https://gtexportal.org/home/datasets). We included 670 samples in our analysis and removed sex chromosomes. The covariates used in eQTL analysis, including top five genotyping principal components (PCs), were obtained from GTEx_Analysis_v8_eQTL_covariates.tar.gz, which was downloaded from GTEx portal (https://gtexportal.org/home/datasets). Then, we further performed a probabilistic estimation of expression residuals (PEER) analysis to adjust for top five genotyping PCs, age, and other potential confounding factors (PEERs)[4] for downstream prediction model building. There is a description of how to download and use the PEER tool here: https://github.com/PMBio/peer/wiki/Tutorial.

**1.1.2	genotype file:**  
The whole genome sequencing file, GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf, was downloaded from dbGaP (https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v8.p2). The genotype dataset is quality controlled using the tool PLINK[5] (https://zzz.bwh.harvard.edu/plink/ ). Multiple QC steps were applied by excluding variants with missingness rate > 0.1, high deviations from Hardy-Weinberg equilibrium at p<10-6, and removing samples with missingness rate > 0.1. To be noticed that, we include all common variants, low frequency variants and rare variants. 

The input genotype file ("genotype_file") with the format as below:

 CHR,LOC,ID-0,ID-2,ID-3,ID-4,ID-5, … … \
 1,933303,0,0,0,0,0,0,0,0,1,0,0,0, … … \
 1,933411,1,2,2,2,2,2,2,2,2,2,2,2, … … \
 1,933653,0,0,0,0,0,0,0,0,0,0,0,0, … … \
 … …

**1.1.3	SNP annotation file:** \
The input snp annotation file ("snp_annot_file"), contains annotations for all variants with the format as below (the same format as bim file): 
1       chr1_10291_C_T_b38      0       10291   T       C \
1       chr1_10291_C_T_b38      0       10291   T       C \
1       chr1_10419_T_C_b38      0       10419   C       T \
1       chr1_10423_C_G_b38      0       10423   G       C \
… …


**1.1.4	Gene annotation file:** \
The input gene annot file ("gene_annot_file") is downloaded from GENCODE: https://www.gencodegenes.org/human/release_26.html, in the GTF format and build in GRCh38.


**1.1.5	GWAS file:** \
The input GWAS file ("gwas_file") contains "CHR" and "LOC" columns, we just need to make sure that all the SNPs being trained in the SuSiE model can be found in the GWAS dataset.

#### 1.2. Training SuSiE model:
We processed one chromosome at a time by executing this code, take chromsome 1 as an example:\
`Rscript ./CODE/Susie_Gene_Chr.R  1`

### 2. Feature aggregation using SKAT:
Running the command:

`java -jar rvTWAS.jar rvTWAS -format csv|plink -input_genotype INPUT_GENOTYPE_FILE -input_phenotype INPUT_PHENOTYPE_FILE -input_phenotype_column INPUT_PHENOTYPE_COLUMN_START:2|6 -input_phenotype_type PHENOTYPE_TYPE: continuous|binary -snp_info_path WEIGHT_FILE  -pheno_id INPUT_GENE_ID  -plink PLINK_BINARY_FILE_PATH  -Rscript RSCRIPT_BINARY_FILE_PATH -output_folder OUTPUT_FOLDER_PATH`

A simple example are described below. Users can get the final p-value result under the folder: OUTPUT_FOLDER_PATH. Take the function "LMM-Kernel" and IDP "25904-2.0" as an example:

#### 2.1. If trying csv format, the command line is:

`java -jar ./CODE/rvTWAS.jar rare-TWAS -format csv -input_genotype ./EXAMPLE/CSV_FORMAT/example.csv -input_phenotype ./EXAMPLE/CSV_FORMAT/example.tsv -input_phenotype_column 2 -input_phenotype_type binary -snp_info_path ./EXAMPLE/Susie_weights.txt -pheno_id ENSG00000239961.2 -plink /PATH/TO/plink -Rscript /PATH/TO/Rscript -output_folder /PATH/TO/OUT_FOLDER`

#### 2.2. If users want to try plink format, the command line is:

`java -jar ./CODE/rvTWAS.jar rare-TWAS -format plink -input_genotype ./EXAMPLE/PLINK_FORMAT/example.tped -input_phenotype ./EXAMPLE/PLINK_FORMAT/example.tfam -input_phenotype_column 6 -input_phenotype_type binary -snp_info_path ./EXAMPLE/Susie_weights.txt -pheno_id ENSG00000239961.2 -plink /PATH/TO/plink -Rscript /PATH/TO/Rscript -output_folder /PATH/TO/OUT_FOLDER`

Please note that, to consistent with plink format, the phenotype is set to missing (normally represented by -9) if unspecified. It must be a numeric value. Case/control phenotypes are normally coded as control = 1, case = 2.The rvTWAS.result under /PATH/TO/OUT_FOLDER is the final output file by rvTWAS.

### References
[1]Wang, G., A. Sarkar, P. Carbonetto and M. Stephens, 2020 A simple new approach to variable selection in regression, with application to genetic fine mapping. Journal of the Royal Statistical Society Series B-Statistical Methodology 82: 1273-1300.\
[2]Wu, M. C., S. Lee, T. Cai, Y. Li, M. Boehnke et al., 2011 Rare-variant association testing for sequencing data with the sequence kernel association test. Am J Hum Genet 89: 82-93.\
[3]Lee, S., M. J. Emond, M. J. Bamshad, K. C. Barnes, M. J. Rieder et al., 2012 Optimal unified approach for rare-variant association testing with application to small-sample case-control whole-exome sequencing studies. Am J Hum Genet 91: 224-237.\
[4]Stegle, O., L. Parts, M. Piipari, J. Winn and R. Durbin, 2012 Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. Nat Protoc 7: 500-507.\
[5]Purcell, S., B. Neale, K. Todd-Brown, L. Thomas, M. A. Ferreira et al., 2007 PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 81: 559-575.

### Contacts
Jingni He: jingni.he1@ucalgary.ca<br>
Qingrun Zhang: qingrun.zhang@ucalgary.ca<br>

### Copyright License (MIT Open Source)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
