## 
##       @file  gemma_on_plink.py
##      @brief  Taking the output of multipleDB allready in plink bed/bim format and cut
##              the full relatedness matrix (kinship) to only the relevant accessions.
##              Also create a fam file of the phenotypes and then run GEMMA on this files
##              to get an association score for each k-mer corrected for population structure.
## 
## Detailed description starts here.
## 
##     @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
## 
##   @internal
##     Created  07/31/18
##    Revision  $Id: doxygen.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
##     Company  Max Planck Institute for Developmental Biology Dep 6
##   Copyright  Copyright (c) 2018, Yoav Voichek
## 
## This source code is released for free distribution under the terms of the
## GNU General Public License as published by the Free Software Foundation.
## =====================================================================================
## 
import os
from glob import glob
gemma_path = "/ebio/abt6/yvoichek/smallproj/prefix/bin/gemma"
def get_file_type_in_dir(dir_name, type_suffix):
    fns = glob("%s/*.%s" % (dir_name, type_suffix))
    if len(fns) != 1:
        raise Exception("File count not = 1 %s %d\n%s"  % (type_suffix, len(fns), dir_name))
    return fns[0]

def get_column(fn, index, sep = "\t"):
    return [line.split(sep)[index] for line in file(fn,"r").read().split("\n")[:-1]]

def create_new_relatedness_matrix(fn_in, fn_out, indices):
    def cut_ind(x, ind):
        return [x[ind[i]] for i in range(len(ind))]
    m = [cut_ind(x.split("\t")[:-1],indices) for x in file(fn_in,"r").read().split("\n")[:-1]]
    m = cut_ind(m, indices)
    # convert back to text
    fout = file(fn_out, "w")
    # in the original file all lines ended with "\t" and there was an empty last line
    # I am trying to keep this standard
    fout.write("\n".join(["\t".join(x)+"\t" for x in m]) + "\n")
    fout.close()

# First version will be set with specific params, but then we will correct to be more general
## General parameters:
base_path = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/"
kinship_matrix = base_path + "tests/build_kinship_matrix_of_all_1001G/kinshipname_all_miss_0_75.cXX.txt"
## to have the order of accession indices
original_fam = base_path + "tests/build_kinship_matrix_of_all_1001G/1001genomes_snp-short-indel_only_ACGTN.vcf.plink.fam.original"


## Specific run parameters:
data_dir = base_path + "tests/run_gemma/43_avrRpt2/"
pheno_fn = get_file_type_in_dir(data_dir, "int_1001g")
bed_fn   = get_file_type_in_dir(data_dir, "bed")
bim_fn   = get_file_type_in_dir(data_dir, "bim")
base_name = bim_fn[:-4]
fam_fn   = base_name + ".fam"

# new file for kinship matrix
kinship_fn = data_dir + "kinship_matrix"
## create a sub matrix of the relevant accessions
orig_acc = get_column(original_fam, 0, " ")
current_acc = get_column(pheno_fn, 0)[1:]
print "We have %d accession in the kinship matrix" % len(orig_acc)
print "We have %d in this phenotype" % len(current_acc)

if len(current_acc) != len([x for x in current_acc if x in orig_acc]):
    raise Exception("Not all the accession in the phenotype also found in the kinship matrix")


index_to_cut = [] ## Finding the indices of current accessions in the general list
for i in range(len(current_acc)):
    index_to_cut.append(orig_acc.index(current_acc[i]))


create_new_relatedness_matrix(kinship_matrix, kinship_fn, index_to_cut)

# building the fam file to run gemma
cur_cmd = r"""cat %s | tail -n +2 | awk '{print $1 " " $1 " 0 0 0 " $2}' > %s""" % (pheno_fn, fam_fn)
print "run cmd: ", cur_cmd 
os.system(cur_cmd)

# run gemma on all k-mers
cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s &" % \
        (gemma_path, base_name,  kinship_fn, data_dir + "output", "kmers")
print "run cmd: ", cur_cmd
os.system(cur_cmd)

# Organizing files to run gemma using the 1001g snps
full_matrix_path = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/tests/build_kinship_matrix_of_all_1001G/"
bed_1135_fn = full_matrix_path + "1001genomes_snp-short-indel_only_ACGTN.vcf.plink.bed"
bim_1135_fn = full_matrix_path + "1001genomes_snp-short-indel_only_ACGTN.vcf.plink.bim"

new_1135 = "1135genomes.plink"
new_bed_1135 = data_dir + new_1135 + ".bed"
new_bim_1135 = data_dir + new_1135 + ".bim"
new_fam_1135 = data_dir + new_1135 + ".fam"

os.system("cp %s %s" % (bed_1135_fn, new_bed_1135))
os.system("cp %s %s" % (bim_1135_fn, new_bim_1135))

#building a fam file to run vs. all snps (so it will be expanded to the 1135 format)
orig_fam_acc = get_column(original_fam, 0, " ")
cur_fam_acc = get_column(fam_fn, 0, " ")
cur_fam_pheno = get_column(fam_fn,5, " ")
fout = file(new_fam_1135, "w")
for (i, acc_n) in enumerate(orig_fam_acc):
    if acc_n in cur_fam_acc:
        fout.write("%s %s 0 0 0 %s\n" % (acc_n, acc_n, cur_fam_pheno[cur_fam_acc.index(acc_n)]))
    else:
        fout.write("%s %s 0 0 0 -9\n" % (acc_n, acc_n))
fout.close()

# run gemma (using the full kinship matrix ofcourse...)
cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s" % \
        (gemma_path, data_dir + new_1135,  kinship_matrix, data_dir + "output", "snp_1135G")
print "run cmd: ", cur_cmd
os.system(cur_cmd)

# should delete coppied files in the end (!)
