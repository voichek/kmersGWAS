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

## General parameters:
import subprocess
import time
base_path = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/"

execfile(base_path + "scripts/general_params.py")

##gemma_path = "/ebio/abt6/yvoichek/smallproj/prefix/bin/gemma-0.98-linux-static"
gemma_path = "/ebio/abt6/yvoichek/smallproj/prefix/bin/gemma"
YV_corr = base_path + "bin/F_correlate_kmers_to_phenotype"
grammar_transformation = base_path + "scripts/R_transform_phenotype/transform_and_permute_phenotypes.R"
kinship_intersect_script = base_path + "scripts/align_kinship_phenotype.py"
DBs_list_fn = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/DBs/1001G_31mers/DBs_to_use"
PERM_N = 100
PARALLEL = 25 
RES_KMERS = 1000001
def count_running_gemma():
    return int(subprocess.check_output("ps -c | grep gemma | wc -l" , \
            shell=True, stderr = subprocess.PIPE)[:-1])

def get_best_pval(fn):
    return subprocess.check_output("cat %s | tail -n +2 | cut -f9 | awk '$1 < 0.001' | sort -g | head -n1"\
            % fn, shell=True, stderr = subprocess.PIPE)[:-1]

def main():
    ## Creating the directory to work in
    name_run = sys.argv[1]
    data_dir = base_path + "tests/pipeline_GrammarGamma/%s/" % name_run
    os.system("mkdir %s" % data_dir)
    print "working on: ", data_dir
    
    # Copying the phenotype data to the directory
    print "Copying phenotype data to directory"
    orig_phenotype_data = glob(base_path + "/tests/phenotype_example/" + name_run + ".pheno.int_1001g")
    if len(orig_phenotype_data) != 1: # has to be zero as glob was specific
        raise Exception("Couldn't find phenotype file")
    os.system("cp %s %s" % (orig_phenotype_data[0], data_dir))

    # Creating kinship matrix
    phenotype_copied_fn = get_file_type_in_dir(data_dir, "int_1001g")
    cur_cmd = "python2.7 %s %s %s %s"
    cur_cmd = cur_cmd % (kinship_intersect_script, phenotype_copied_fn, data_dir + "/" + name_run, DBs_list_fn)
    os.system(cur_cmd)

    # Transformation of the phenotypes and creating permutations of phenotype
    kinship_fn_emma = get_file_type_in_dir(data_dir, "kinship.emma")
    kinship_fn_gemma = get_file_type_in_dir(data_dir, "kinship.gemma")
    pheno_fn = get_file_type_in_dir(data_dir, "pheno")

    cur_cmd = "Rscript %s %s %s %d %s"
    cur_cmd = cur_cmd % (grammar_transformation, pheno_fn, kinship_fn_emma, PERM_N, data_dir + "/" + name_run)
    os.system(cur_cmd)

    print "run correlation with k-mers"
    perm_pheno_fn = get_file_type_in_dir(data_dir, "phenotype")
    trans_pheno_fn = get_file_type_in_dir(data_dir, "phenotype.trans")
    output_F_corr = "%s/F_corr_output" % data_dir
    os.system("mkdir %s" % output_F_corr)

    cur_cmd = "%s -p %s -b %s -o %s -n %d --parallel %d --paths_file %s --kmers_file %s"  %\
            (YV_corr, trans_pheno_fn, name_run, output_F_corr, RES_KMERS, PARALLEL, DBs_list_fn, "kmers_31mers_5min_1008acc_sorted")
    print cur_cmd
    os.system(cur_cmd)

    # for all the fam files, move to fam.orig and create a new fam file with the non-transformed phenotype
    phenotypes_names = file(perm_pheno_fn,"r").read().split("\n")[0].split("\t")[1:]
    print "We have %d phenotypes" % len(phenotypes_names)
    for (p_ind, p_name) in enumerate(phenotypes_names):
        bed_fn = get_file_type_in_dir(output_F_corr, "%s.bed" % p_name)
        bim_fn = get_file_type_in_dir(output_F_corr, "%s.bim" % p_name)
        fam_fn = get_file_type_in_dir(output_F_corr, "%s.fam" % p_name)
        os.system("mv %s %s.orig" % (fam_fn, fam_fn))
        base_name = bed_fn[:-4]
        # building the fam file to run gemma
        cur_cmd = r"""cat %s | tail -n +2 | awk '{print $1 " " $1 " 0 0 0 " $%d}' > %s""" % \
                (perm_pheno_fn, p_ind +2,  fam_fn)
        print cur_cmd 
        os.system(cur_cmd)
        # run gemma on best k-mers
        cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -maf 0.05 -miss 0.5 &" % \
                (gemma_path, base_name,  kinship_fn_gemma, output_F_corr + "/output", p_name)
        while count_running_gemma() >= (PARALLEL-1):
            time.sleep(30)
        
        time.sleep(1)
        print "run cmd: ", cur_cmd
        os.system(cur_cmd)

    # Organizing files to run gemma using the 1001g snps
    full_matrix_path = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/tests/build_kinship_matrix_of_all_1001G/"
    bed_1135_fn = full_matrix_path + "1001genomes_snp-short-indel_only_ACGTN.vcf.plink.bed"
    bim_1135_fn = full_matrix_path + "1001genomes_snp-short-indel_only_ACGTN.vcf.plink.bim"
    
    dir_1001G_vcf = data_dir + "1001G"
    os.system("mkdir " +  dir_1001G_vcf)
    new_1135 = "/1135genomes.plink"
    new_bed_1135 = dir_1001G_vcf + new_1135 + ".bed"
    new_bim_1135 = dir_1001G_vcf + new_1135 + ".bim"
    new_fam_1135 = dir_1001G_vcf + new_1135 + ".fam"
####    
    os.system("cp %s %s" % (bed_1135_fn, new_bed_1135))
    os.system("cp %s %s" % (bim_1135_fn, new_bim_1135))
    
    #building a fam file to run vs. all snps (so it will be expanded to the 1135 format)
    orig_fam_acc = get_column(original_fam, 0, " ")
    cur_fam_acc = get_column(perm_pheno_fn, 0, "\t")[1:]
        
    cur_fam_pheno = [get_column(perm_pheno_fn,i+1, "\t")[1:] for i in range(PERM_N+1)]
    fout = file(new_fam_1135, "w")
    for (i, acc_n) in enumerate(orig_fam_acc):
        fout.write("%s %s 0 0 0" % (acc_n, acc_n))
        if acc_n in cur_fam_acc:
            for p_ind in range(PERM_N+1):
                fout.write(" %s" % cur_fam_pheno[p_ind][cur_fam_acc.index(acc_n)])
        else:
            for p_ind in range(PERM_N+1):
                fout.write(" -9")
        fout.write("\n")
    fout.close()
    # run gemma (using the full kinship matrix ofcourse...)
    for (p_ind, p_name) in enumerate(phenotypes_names):
        cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -n %d -maf 0.05 -miss 0.5 &" % \
                (gemma_path, dir_1001G_vcf+ new_1135,  kinship_matrix_GEMMA, dir_1001G_vcf + "/output", \
                p_name, p_ind+1)
        while count_running_gemma() >= (PARALLEL-1):
            time.sleep(30)
        
        time.sleep(1)
        print "run cmd: ", cur_cmd
        os.system(cur_cmd)
    
    while count_running_gemma() != 0:
        time.sleep(30)
    best_pvalues_fn = data_dir + "best_pvals.csv"
    fout = file(best_pvalues_fn,"w")

    # condense all the best p_values to one file
    output_kmers = output_F_corr + "/output"
    output_snps = dir_1001G_vcf + "/output"
    fout.write("name\tk_mers\tSNPS\n")
    for (p_ind, p_name) in enumerate(phenotypes_names):
        fn_kmers = output_kmers + "/%s.assoc.txt" % p_name
        fn_snps  = output_snps + "/%s.assoc.txt" % p_name
        fout.write("%s\t%s\t%s\n" % (p_name, get_best_pval(fn_kmers),get_best_pval(fn_snps)))

    fout.close()



if __name__ == "__main__":
    main()
