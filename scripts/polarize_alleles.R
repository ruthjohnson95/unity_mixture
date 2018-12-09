#!/usr/bin/env Rscript

library(stringr)
library(data.table)
library(dplyr)
library("optparse")
 
#option_list = list(
#	      make_option(c("--gwas_file"), type="character", metavar="character"),
#	      make_option(c("--mrca_chr_file"), type="character", metavar="character"),
#	      make_option(c("--chr"), type="int", metavar="int"),
#	      make_option(c("--outdir"), type="character", metavar="character")
#); 
 
option_list = list(
  make_option(c("--gwas_file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
	      make_option(c("--mrca_chr_file"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
	      make_option(c("--outdir"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
              make_option(c("--chr"), type="integer",
              help="output file name [default= %default]", metavar="character")
); 
 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

numextract <- function(string){ 
  as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
} 

# gwas data
gwas_file=opt$gwas_file
mrca_file=opt$mrca_chr_file # mrca file by chr 
chr=opt$chr
outdir=opt$outdir 

gwas=fread(gwas_file, header=TRUE) # assumes SNP, CHR, A1

mrca=fread(mrca_file, header=F)
colnames(mrca) <- c("CHR", "BP", "AA")

# only use chromosome number 
mrca$CHR <- numextract(mrca$CHR)

# only look at gwas chri
chr_inds = which(gwas$CHR == chr)
gwas_i = gwas[chr_inds]

# convert to uppercase bc fread is switches alleles to lowercase 
mrca$AA <- lapply(mrca$AA,toupper)
known_AA_inds<-which(mrca$AA != '-')
mrca <- mrca[known_AA_inds, ]

# use all common variables across two tables (SNP, CHR)
gwas_i = inner_join(gwas_i, mrca, by=NULL)

# filter out SNPs with unknown AA 
##known_AA_inds <- which(gwas_i$AA != '-')
#gwas_i = gwas_i[known_AA_inds]
	    
# loop through and flip alleles 
gwas_i$BETA_STD<- ifelse(gwas_i$AA == gwas_i$A1, gwas_i$BETA_STD, gwas_i$BETA_STD*(-1))

# save gwas
trait=basename(sub('\\.txt$', '', gwas_file) )
chr_str=paste('chr', chr, sep='')
outfile=file.path(trait, chr_str, 'mrca.txt',fsep='_')
#write.table(gwas_i, file.path(outdir, outfile), quote=F, row.names=F)
final_out_file=file.path(outdir, outfile)
fwrite(gwas_i, final_out_file, sep=' ')

