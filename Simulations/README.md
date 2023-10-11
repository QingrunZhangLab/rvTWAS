## Overview 
Aligning to the rvTWAS protocol borrowing gene expression from a reference dataset, the simulations generated two datasets: the gene expressions in the reference genotype-expression dataset and the phenotype in the main GWAS dataset (Figure 1C). Note that in both datasets, genotypes are not simulated. 
Based on the GTEx genotype, we simulated gene expressions in the reference dataset (Figure 1C, upper) corresponding to the analysis in Step 1 of any TWAS protocol. Based on the 1000 Genomes Project genotype, we simulated phenotype in the GWAS dataset (Figure 1C, lower) corresponding to the analysis in Step 2 of any TWAS analysis. 

### Step1: 
1) Following the conventional heterogeneity model, in which any rare variant contributes to the phenotypic changes independently, we form “pseudo-SNVs” to represent the aggregated effects of multiple rare variants. More specifically, a number of (M = 0, 50, 100, or 200) rare variants (with MAF <0.5%) in the focal gene region were randomly selected and combined into a pseudo-SNV. The pseudo-SNV is coded by aggregation: the subjects carrying at least one of the M rare variants were coded as 1; and the subjects that do not carry any rare variant were coded as 0. 
2) Following an additive model, we counted the number of rare variants carried by an individual as the contribution of the rare variants to its phenotype, which is equivalent to a “Sum” of these rare variants. 

Using target gene id = ENSG00000239961.2 as an example:

`python ./CODES/WriteInfoFile_Step1.py ENSG00000239961.2`

### Step2:
Acknowledging that common variant may also contribute to gene expression or phenotype, the overall genetic component is modelled as the sum of this rare variant effect (in the form of either pseudo-SNV or Sum-effect) and common variants (Figure 1D lower). More specifically, five common variants (with MAF >5%) were randomly selected as causal variants, and the weighted sum of these six causal variants (5 common + 1 pseudo/sum SNVs) will be considered as the genetic component (of the expression or phenotype). The weight of each causal variant is sampled from a standard normal distribution N(0,1). 

Using target gene id = ENSG00000239961.2 and M=200 as an example:

For “pseudo-SNVs”:

`python ./CODES/WriteInfoFile_Step2.py ENSG00000239961.2 200`

For “Sum”:

`python ./CODES/WriteInfoFile_Step2_SUM.py ENSG00000239961.2 200`

### Step3:
Two typical genetic architectures representing relationships between genotype, expression and phenotype are considered: causality, where genotypes alter phenotype via the expression (Figure 1E upper), and pleiotropy, where genotypes contribute to phenotype and expression independently (Figure 1E lower). Under both scenarios, we simulated expression and phenotypes using an additive genetic model, in which phenotypes and expression are caused by a weighted sum of genetic effects. In all cases, we simulated a genetic component first and then incorporate it to the expression or phenotypes using prespecified values of genetic component, i.e., heritability.

Given the case under “pseudo-SNVs” model, scenario=causality, target gene id = ENSG00000239961.2, M=200, phenotype heritability=0.1, and expression heritability=0.15, as an example:

Using GTEx genotype:

`java -jar ./CODES/eQTL.jar sim-phenotype ENSG00000239961.2_gtex_genotype_pseudorare_200.csv INFO_FUNC_200_ENSG00000239961.2.txt causality 0.15 0.1 EXP_Rare200_Ex0.15_Ph0.1_causality_ENSG00000239961.2.txt PHENO_Rare200_Ex0.15_Ph0.1_causality_ENSG00000239961.2.txt`

Using1000 Genomes Project genotype:

`java -jar ./CODES/eQTL.jar sim-phenotype ENSG00000239961.2_hg1000_genotype_pseudorare_200.csv INFO_FUNC_200_ENSG00000239961.2.txt causality 0.15 0.1 EXP_Rare200_Ex0.15_Ph0.1_causality_ENSG00000239961.2.txt PHENO_Rare200_Ex0.15_Ph0.1_causality_ENSG00000239961.2.txt`
