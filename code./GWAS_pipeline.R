#GWAS PIPELINE
library(optparse)
library(data.table)
library(dplyr)
set.seed(10)
option_list <- list(
  make_option(c("-v", "--vcf_gz"), type="character", default=FALSE,
    help="Insert path to vcf.gz file"),
  make_option(c("-o", "--out_vcf"), type="character", default=FALSE,
    help="Name of filtered vcf.gz to output "),
  make_option(c("-f", "--filter_vcf"), type="character", default=TRUE,
    help="Filter vcf file: arguments should be TRUE or FALSE"),
  make_option(c("-g", "--make_geno"), type="character", default=FALSE,
    help="Make genotype file: arguments should be TRUE or FALSE"),
  make_option(c("-e", "--out_geno"), type="character", default=TRUE,
    help="Name of geno file output"), 
  make_option(c("-n", "--sample_names"), type="character", default=TRUE,
    help="Text file with sample names"),    
  make_option(c("-c", "--covar_file"), type="character", default=TRUE,
    help="Name of covariate file"),
  make_option(c("-k", "--k_file"), type="character", default=TRUE,
    help="Name of file containing kinship matrix"),
  make_option(c("-l", "--lineages"), type="character", default=TRUE,
    help="If lineages==all, all samples will be used"),
  make_option(c("-d", "--drug"), type="character", default=TRUE,
    help="Name of drug (phenotypes)"),
  make_option(c("-p", "--pheno_file"), type="character", default=TRUE,
    help="Name of phenotype file to use"),
  make_option(c("-r", "--results_file"), type="character", default=TRUE,
    help="Name of results file to output"),
  make_option(c("-x", "--perform_gwas"), type="character", default=FALSE,
    help="TRUE or FALSE"),
  make_option(c("-z", "--make_pheno"), type="character", default=FALSE,
    help="TRUE or FALSE"),
    make_option(c("-y", "--metadata_file"), type="character", default=FALSE,
    help="Name of metadata file")
    )
#parse options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
#Functions
#Function to filter vcf by gene region
filter_vcf <- function(vcf_gz, out_vcf){
  cat(" *Filtering vcf using region using the following files \n")
  #print(region)
  #print(vcf_gz)
  #system(paste("bcftools view -Oz -R", region, vcf_gz, ">" , out_vcf))
  system(paste("bcftools norm -Oz -m-any", vcf_gz, ">", out_vcf))}
#Function to make and format genotype file for further analyses
make_geno <- function(out_vcf, geno,  sample_names){
  cat(" *Making geno file using the following files \n")
  system(paste("bcftools index -c --force", out_vcf))
  system(paste("bcftools norm -m -any -Oz", out_vcf, "> geno.vcf.gz"))
  system(paste("bcftools query -f '%CHROM\t%POS\t%ALT[\t%GT]\n' geno.vcf.gz | tr '|' '/' | sed 's/\\.\\/\\./Ns/g' | sed 's/0\\/1/0.5/g' | sed 's/[123456789]\\/[123456789]/1/g' | sed 's/0\\/0/0/g' >", geno))
  system(paste("bcftools query -l geno.vcf.gz >", sample_names))
  system(paste("awk \'BEGIN{OFS=\"_\"} {print $1,$2,$3}\'", geno, "> tmp"))
  system(paste("awk \'NR==FNR{a[NR]=$0;next}{$1=a[FNR]}1\' tmp ", geno, ">tmp2 &&  mv tmp2", geno, "&& rm tmp2 && rm tmp"))  
  system(paste("sed -i 's/ /, /g'", geno))
  system(paste("sed -i 's/N/NA/g'", geno))
  system(paste("gzip", geno))
   }
#Function to make pheno file
make_pheno <- function(metadata_file, drug, lineages, pheno_file){
if (lineages=="all"){
metadata <- read.csv(metadata_file, header=TRUE)
pheno <- metadata[,drug]
write.table(pheno, pheno_file, col.names=FALSE, row.names=FALSE, quote = FALSE)
}
else {
metadata <- read.csv(metadata_file, header=TRUE)
metadata$lineage <- gsub("lineage1.*", "lineage1", metadata$lineage)
metadata$lineage <- gsub("lineage2.*", "lineage2", metadata$lineage)
metadata$lineage <- gsub("lineage3.*", "lineage3", metadata$lineage)
metadata$lineage <- gsub("lineage4.*", "lineage4", metadata$lineage)
metadata$lineage <- gsub("lineage5.*", "lineage5", metadata$lineage)
metadata$lineage <- gsub("lineage6.*", "lineage6", metadata$lineage)
metadata$lineage <- gsub("lineage7.*", "lineage7", metadata$lineage)
if(metadata[,"lineage"]!=lineages){
metadata[,drug] <- "NA"
}
pheno <- metadata[,drug]
write.table(pheno, pheno_file, col.names=FALSE, row.names=FALSE, quote = FALSE)
}
}
#Function to perform GWAS using Gemma
perform_gwas <- function(geno, pheno_file, covar_file, k_file, outfile){
cat(" * Compute p-vals. \n")
geno_gz <- paste0(geno,".gz")
cat(geno_gz)
system(paste("gemma -g", geno_gz," -p", pheno_file, "-k", k_file, "-c", covar_file,  "-lmm 1 -maf 0.001 -o", outfile),  ignore.stdout = FALSE)
}
#Run
if(opt$filter_vcf==TRUE){
filter_vcf(opt$vcf_gz, opt$out_vcf)
}
if(opt$make_geno==TRUE){
make_geno(opt$out_vcf, opt$out_geno,  opt$sample_names)
}
if(opt$make_pheno==TRUE){
make_pheno(opt$metadata_file, opt$drug, opt$lineages, opt$pheno_file)
}
if(opt$perform_gwas==TRUE){
perform_gwas(opt$out_geno, opt$pheno_file, opt$covar_file, opt$k_file, opt$results_file)
}


