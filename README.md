# GWAS
Genome wide association analysis to identify novel drug-resistance mutations and loci from M. tuberculosis WGS
Code is supplied in pipelines which can be run using the following examples. 
There are a few requirements:
* A vcf file containing the genotype data for samples.
* A metadata file. The first column must contain the sample id with the name 'id', there should also be the column 'lineage' and name of drug. 
* Covariates file (tsv)- contains covariates used in GWAS analysis.

## Examples of how to run code on command line
### GWAS
```
Rscript GWAS_pipeline.R --vcf_gz "tb.vcf.gz" --out_vcf "gwas_snp.vcf.gz" --filter_vcf FALSE --sample_names "sample_names.txt" --make_geno FALSE --out_geno "tb.geno"  --covar_file "covar.tsv" --k_file "output/tb.cXX.txt" --lineages "all" --drug "rifampicin" --pheno_file "rif_gwas_phen.txt" --results_file "rifampicin.txt" --perform_gwas TRUE --make_pheno TRUE --metadata_file "metadata.csv"
```
### Permutations
```
Rscript Permutations.R --pheno_name "rifampicin" --geno_file "tb.geno" --pheno_file "rif_gwas_phen.txt" --covars "covar.tsv"
```
