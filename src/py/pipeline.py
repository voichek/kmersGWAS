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
import sys

paths = {}
paths["base_path"] = "/".join(sys.path[0].split("/")[:-2]) + "/"
paths["gen_script"] = paths["base_path"] + "src/py/functions.py"
paths["assoicte_kmers"] = paths["base_path"] + "bin/F_correlate_kmers_to_phenotype"
paths["Permute_phenotype"] = paths["base_path"] + "src/R/transform_and_permute_phenotypes.R"
paths["kinship_intersect_script"] = paths["base_path"] + "src/py/align_kinship_phenotype.py"
paths["parser_script"] = paths["base_path"] + "src/py/pipeline_parser.py"
paths["average_pheno_script"] = paths["base_path"] + "src/awk/average_phenotypes.awk"

####################################################################################################
execfile(paths["gen_script"])
execfile(paths["parser_script"])
####################################################################################################

def main():
    args = parser.parse_args()

    ## Creating the directory to work in
    print args
    print args.outdir
    if dir_exist(args.outdir):
        raise Exception("Output directory allready exists")

    os.system("mkdir %s" % args.outdir)
    
    # Copying the phenotype data to the directory
    if args.verbose:
        print "Copying phenotype data to directory"
    
    paths["pheno_orig_fn"] = "%s/%s.original_pheno" % (args.outdir, args.name)
   # os.system("cp %s %s" % (get_file(args.fn_phenotype), paths["pheno_orig_fn"]))
    os.system("cat %s | tail -n +2 | awk -f %s > %s" % \
            (get_file(args.fn_phenotype),paths["average_pheno_script"],paths["pheno_orig_fn"]))


    # Information on SNPs dataset
    paths["snps_fam"] = args.snps_matrix + ".fam"


    # Creating kinship matrix
    paths["pheno_intersected_fn"] = "%s/%s.intersected_phenoo" % (args.outdir, args.name)
    paths["kinship_ibs_fn"] = "%s/%s.kinship.ibs" % (args.outdir, args.name)
    paths["kinship_gemma_fn"] = "%s/%s.kinship.gemma" % (args.outdir, args.name)
    
    cur_cmd = "python2.7 %s --pheno %s --fam_file %s --kinship_file %s --output_pheno  %s " + \
            "--output_kinship %s --DBs_list %s"
    cur_cmd = cur_cmd % (paths["kinship_intersect_script"], paths["pheno_orig_fn"], paths["snps_fam"], 
            args.kinship_IBS, paths["pheno_intersected_fn"], paths["kinship_ibs_fn"], args.db_list)
    if args.verbose:
        cur_cmd = cur_cmd + " -v"
        print "RUN: ", cur_cmd

    os.system(cur_cmd)
    
    cur_cmd = "python2.7 %s --pheno %s --fam_file %s --kinship_file %s --output_kinship %s --DBs_list %s"
    cur_cmd = cur_cmd % (paths["kinship_intersect_script"], paths["pheno_orig_fn"], paths["snps_fam"], \
            args.kinship_GEMMA, paths["kinship_gemma_fn"], args.db_list)
    if args.verbose:
        cur_cmd = cur_cmd + " -v"
        print "RUN: ", cur_cmd
    os.system(cur_cmd)

    # Transformation of the phenotypes and creating permutations of phenotype
    paths["pheno_permuted_fn"] = "%s/%s.permuted_pheno" % (args.outdir, args.name)
    paths["pheno_permuted_transformed_fn"] = "%s/%s.permuted_transformed_pheno" % (args.outdir, args.name)
    paths["EMMA_perm_log_fn"] = "%s/EMMA_perm.log" % (args.outdir)

    cur_cmd = "Rscript %s %s %s %d %s %s %s"
    cur_cmd = cur_cmd % (paths["Permute_phenotype"] , paths["pheno_intersected_fn"], paths["kinship_ibs_fn"],\
           args.n_permutations, paths["pheno_permuted_fn"], paths["pheno_permuted_transformed_fn"], paths["EMMA_perm_log_fn"])
    if args.verbose:
        print "RUN: ", cur_cmd
    os.system(cur_cmd)

    paths["kmers_associations_dir"] = "%s/F_corr_output" % args.outdir
    os.system("mkdir %s" % paths["kmers_associations_dir"])

    cur_cmd = "%s -p %s -b %s -o %s -n %d --parallel %d --paths_file %s --kmers_table %s --kmer_len %d"  %\
            (paths["assoicte_kmers"], paths["pheno_permuted_transformed_fn"], args.name, 
                    paths["kmers_associations_dir"], args.n_kmers, args.parallel, args.db_list, 
                    args.kmers_table, args.kmers_len)
    if args.verbose:
        print cur_cmd
    os.system(cur_cmd)

    # for all the fam files, move to fam.orig and create a new fam file with the non-transformed phenotype
    phenotypes_names = file(paths["pheno_permuted_fn"] ,"r").read().split("\n")[0].split("\t")[1:]
    if args.verbose:
        print "We have %d phenotypes" % len(phenotypes_names)
    for (p_ind, p_name) in enumerate(phenotypes_names):
        bed_fn = get_file_type_in_dir(paths["kmers_associations_dir"], "%s.bed" % p_name)
        bim_fn = get_file_type_in_dir(paths["kmers_associations_dir"], "%s.bim" % p_name)
        fam_fn = get_file_type_in_dir(paths["kmers_associations_dir"], "%s.fam" % p_name)
        os.system("mv %s %s.orig" % (fam_fn, fam_fn)) # Saving the original fam (which has transformed values)
        base_name = bed_fn[:-4]
        # building the fam file to run gemma
        cur_cmd = r"""cat %s | tail -n +2 | awk '{print $1 " " $1 " 0 0 0 " $%d}' > %s""" % \
                (paths["pheno_permuted_fn"], p_ind +2,  fam_fn)
        if args.verbose:
            print cur_cmd 
        os.system(cur_cmd)
        # run gemma on best k-mers
        cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -maf 0.05 -miss 0.5 &" % \
                (args.gemma_path, base_name,  paths["kinship_gemma_fn"], \
                paths["kmers_associations_dir"] + "/output", p_name)
        while count_running_gemma() >= (args.parallel-1):
            time.sleep(30)
        
        time.sleep(1)
        if args.verbose:
            print "run cmd: ", cur_cmd
        os.system(cur_cmd)

    if args.run_snps:
        paths["snps_associations_dir"] = "%s/SNPs" % args.outdir
        os.system("mkdir %s" % paths["snps_associations_dir"])

        paths["snps_table_fn"] = "%s/snps.plink" % paths["snps_associations_dir"]
        
        os.system("ln -s %s.bed %s.bed" % (args.snps_matrix, paths["snps_table_fn"]))
        os.system("ln -s %s.bim %s.bim" % (args.snps_matrix, paths["snps_table_fn"]))
    
        #building a fam file to run vs. all snps (so it will be expanded to the 1135 format)
        orig_fam_acc = get_column(paths["snps_fam"], 0, " ")
        cur_fam_acc = get_column(paths["pheno_permuted_fn"], 0, "\t")[1:]
        cur_fam_pheno = [get_column(paths["pheno_permuted_fn"],i+1, "\t")[1:] for i in range(args.n_permutations+1)]
        fout = file("%s.fam" % paths["snps_table_fn"], "w")
        for (i, acc_n) in enumerate(orig_fam_acc):
            fout.write("%s %s 0 0 0" % (acc_n, acc_n))
            if acc_n in cur_fam_acc:
                for p_ind in range(args.n_permutations+1):
                    fout.write(" %s" % cur_fam_pheno[p_ind][cur_fam_acc.index(acc_n)])
            else:
                for p_ind in range(args.n_permutations+1):
                    fout.write(" -9")
            fout.write("\n")
        fout.close()
        # run gemma (using the full kinship matrix ofcourse...)
        for (p_ind, p_name) in enumerate(phenotypes_names):
            cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -n %d -maf 0.05 -miss 0.5 &" % \
                    (args.gemma_path, paths["snps_table_fn"], args.kinship_GEMMA, \
                    paths["snps_associations_dir"] + "/output", \
                    p_name, p_ind+1)
            while count_running_gemma() >= (args.parallel-1):
                time.sleep(30)
            time.sleep(1)
            if args.verbose:
                print "run cmd: ", cur_cmd
            os.system(cur_cmd)
        while count_running_gemma() != 0:
            time.sleep(30)
#        best_pvalues_fn = data_dir + "best_pvals.csv"
#        fout = file(best_pvalues_fn,"w")
#
#    # condense all the best p_values to one file
#    output_kmers = output_F_corr + "/output"
#    output_snps = dir_1001G_vcf + "/output"
#    fout.write("name\tk_mers\tSNPS\n")
#    for (p_ind, p_name) in enumerate(phenotypes_names):
#        fn_kmers = output_kmers + "/%s.assoc.txt" % p_name
#        fn_snps  = output_snps + "/%s.assoc.txt" % p_name
#        fout.write("%s\t%s\t%s\n" % (p_name, get_best_pval(fn_kmers),get_best_pval(fn_snps)))
#
#    fout.close()
#
#

if __name__ == "__main__":
    main()
