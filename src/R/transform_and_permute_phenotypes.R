# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# TRANSFORMATION OF PHENOTYPES FOR GRAMMAR-GEMMA SCORE
# Before scoring k-mers we will use the Grammar-Gamma method to transform the phenotypes.
# "normalizing phenotypes" by the relatedeness structure (and we can also include covariates here in the future)
# This program will also have two other functinalities:
# 1. Permutating the phenotypes: to get a significant threshold, the k-mers scoring will be done also on
# permuted phenotype data. As transformation has to be done after permutation, the permutation will also be part
# of this initial pre-proccessing
# 2. To increase speed of the k-mer scoring, scores will be linearly transformed to be positive integers.
# The relative scoring, thus the order of k-mers, should be invariant to linear transfomrations.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# The transformation of the phenotype will be done using code extracted from the GenABEL package. This package is
# not supported anymore - we will take out the code so we can support it in the future.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Regarding relatedness matrix - I am not sure if it matters how I create it from the SNPs file. I tried to 
# create it with the rvtests which also have a Grammar-Gamma implementation, and I got results similiar to the 
# GEMMA kinship matrix, but not exactly the same, and also there was a factor 10 between the two. For now, as not
# to support too much code, and as I would use GEMMA in the next stage for correct p-values, I will use the kinship
# matrix created by GEMMA.
# But to keep the units I will divide by 10 to have the same scaling
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Set environment:
rm(list = ls())
library(MASS)
library(methods)
library(mvnpermute)
library(matrixcalc)

source('/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/src/R/genabel/polygenic.R')
source('/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/src/R/genabel/polylik.R')
source('/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/src/R/emma.R')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
linear_trans_to_natural_numbers <- function(x) {
  x <- x-min(x)
  max_x <- max(x)
  if(max_x > (10^6)) {
    stop('Not clear how to convert phenotype to natural numbers')
  }
  x <- round(x * (10^8 / max_x))
  return(x)
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
args <- commandArgs(trailingOnly = TRUE)

set.seed(123456789) # To have a reporducible results
# Load user parameters
fn_phenotypes <- args[1]
fn_kinship <- args[2]
n_permute <- as.numeric(args[3])
fn_out_phenotypes <- args[4]
fn_out_trans_phenotypes <- args[5]
f_log <- file(args[6])
f_log_grammar <- paste(args[5],".log",sep="")
# 1. Load phenotype file
phenotypes <- read.csv(fn_phenotypes,sep='\t',header=TRUE)
av_pheno <- mean(phenotypes$phenotype_value)
phenotypes$phenotype_value <- phenotypes$phenotype_value-av_pheno
n_acc = nrow(phenotypes)

# 2. Load kinship matrix
K <- read.csv(fn_kinship, sep='\t',header =FALSE)
K <- as.matrix(K[,1:(dim(K)[2]-1)])

### Check matrix is positive semi-definite
if(!is.positive.semi.definite(K)) {
  writeLines('Kinship matrix is not positive semi-definite')
  quit()
}
### 
null <- emma.REMLE(phenotypes$phenotype_value,as.matrix(x = rep(1, n_acc),dim = c(n_acc,1)),K)
herit <- null$vg/(null$vg+null$ve)
COV_MATRIX <- null$vg*K+null$ve*diag(dim(K)[1])
# Logging
writeLines(c(paste('EMMA_n_permutation','=',n_permute), 
             paste('EMMA_n_accessions','=',n_acc), 
             paste('EMMA_vg','=',null$vg), 
             paste('EMMA_ve','=',null$ve), 
             paste('EMMA_herit','=',herit)),f_log)
close(f_log)
if(n_permute > 0) {
  # K_stand <- (n_acc-1)/sum((diag(n_acc)-matrix(1,n_acc,n_acc)/n_acc)*K)*K #Doesn't effect COV_MATRIX

  # Part of second option to permute phenotypes
  # M <- solve(chol(COV_MATRIX))
  # Y_t <- crossprod(M,phenotypes$phenotype_value)
  # M_inv <- ginv(M)
  
  permute_phenotype <- mvnpermute(phenotypes$phenotype_value, rep(1, n_acc), COV_MATRIX, nr=n_permute, seed=123456789)
  
  # 3. Permute phenotypes
  for(i in 1:n_permute) {
    phenotypes[paste("P",i, sep = '')] <- permute_phenotype[,i]
    # cur_perm = sample(1:n_acc)
    # Part of second method to permute phenotypes
    # phenotypes[paste("P",i, sep = '')] <- crossprod(M_inv, Y_t[cur_perm])
  }
}

# 4. Transform phenotypes
trans_phenotypes <- phenotypes
grammar_run_info <- data.frame(phenotype = I(array(data="",dim=n_permute+1)),
                               GRAMMAR_Gamma_Test = array(data="",dim=n_permute+1),
                               GRAMMAR_Gamma_Beta = array(data="",dim=n_permute+1),
                               GRAMMAR_Gamma_esth2 = array(data="",dim=n_permute+1))

for(i in 2:(n_permute+2)) {
  cur_pre_GRAMMAR_Gamma <- polygenic(phenotypes[,i], K,data.frame(c()), quiet = TRUE)
  # Save GRAMMAR Gamma run information
  grammar_run_info$phenotype[i-1] <- colnames(phenotypes)[i]
  grammar_run_info$GRAMMAR_Gamma_Test[i-1] <- cur_pre_GRAMMAR_Gamma$grammarGamma$Test
  grammar_run_info$GRAMMAR_Gamma_Beta[i-1] <- cur_pre_GRAMMAR_Gamma$grammarGamma$Beta
  grammar_run_info$GRAMMAR_Gamma_esth2[i-1] <- cur_pre_GRAMMAR_Gamma$esth2
  trans_phenotypes[,i] <- cur_pre_GRAMMAR_Gamma$pgresidualY
}

# Log GRAMMAR-Gamma run information
write.csv(x = grammar_run_info, file = f_log_grammar, quote=FALSE,sep="\t",row.names = FALSE)
# 5. Linearly transform transformed phenotypes
# trans_phenotypes[,2:(n_permute+2)] <- linear_trans_to_natural_numbers(trans_phenotypes[,2:(n_permute+2)])

# 6. Output phenotypes before and after transformation (4+)5
write.table(x=phenotypes, file=fn_out_phenotypes, eol = "\n",sep = "\t", quote=F, row.names = FALSE)
# options(scipen=20)
write.table(x=trans_phenotypes, file=fn_out_trans_phenotypes, eol = "\n",sep = "\t", quote=F, row.names = FALSE)
# options(scipen=0)  # restore the default
