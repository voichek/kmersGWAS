from glob import glob
import os
import sys

kinship_matrix_GEMMA = base_path + "tests/build_kinship_matrix_of_all_1001G/kinshipname_all_miss_0_75_gk2.sXX.txt"
kinship_matrix_EMMA = base_path + "tests/build_kinship_matrix_of_all_1001G/kinship_emma.maf0_05"
## to have the order of accession indices
original_fam = base_path + "tests/build_kinship_matrix_of_all_1001G/1001genomes_snp-short-indel_only_ACGTN.vcf.plink.fam.original"

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
