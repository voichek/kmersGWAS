## 
##       @file  pipeline.py
##      @brief  A pipeline to associate phenotype with k-mers and snps
## 
## This script contains the logic of associating a phenotype with k-mers or SNPs
## and all the organization and follow of parameters between the different programs and
## scripts used.
## 
##     @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
## 
##   @internal
##     Created  07/31/18
##     Company  Max Planck Institute for Developmental Biology Dep 6
##   Copyright  Copyright (c) 2018, Yoav Voichek
## 
## This source code is released for free distribution under the terms of the
## GNU General Public License as published by the Free Software Foundation.
## =====================================================================================

## General parameters:
import subprocess
import time
import sys

## Save all the relevant directories and files in one dictionary
paths = {}
# scripts/ programs prefix
paths["base_path"] = "/".join(sys.path[0].split("/")[:-2]) + "/"
# script with python functionalities
paths["gen_script"] = paths["base_path"] + "src/py/functions.py"
# program the associate kmers
paths["assoicte_kmers"] = paths["base_path"] + "bin/associate_kmers_with_phenotypes"
# program to associate snps
paths["associate_snps"] = paths["base_path"] + "bin/associate_snps_with_phenotypes"
# script to permute and transform phenotypes
paths["Permute_phenotype"] = paths["base_path"] + "src/R/transform_and_permute_phenotypes.R"
# script to intersect kinship matrix to the current samples
paths["kinship_intersect_script"] = paths["base_path"] + "src/py/align_kinship_phenotype.py"
# logic of parsing the user given params
paths["parser_script"] = paths["base_path"] + "src/py/pipeline_parser.py"
# average phenotypes if there are a few repeats of some accessions
paths["average_pheno_script"] = paths["base_path"] + "src/awk/average_phenotypes.awk"

####################################################################################################
execfile(paths["gen_script"])
execfile(paths["parser_script"])
####################################################################################################


def main():
    ## Read and parse user defined parameters
    args = parser.parse_args()

    ## Creating the directory to work in
    create_dir(args.outdir)

    ## open log file
    paths["log_file"] = "%s/log_file" % args.outdir
    f_log = file(paths["log_file"],"w",buffering=1)
    f_log.writelines(str(args))

    ## Copying the phenotype data to the directory
    paths["pheno_orig_fn"] = "%s/%s.original_phenotypes" % (args.outdir, args.name)
    copy_phenotypes(get_file(args.fn_phenotype), paths["pheno_orig_fn"], paths["average_pheno_script"], f_log)

    ## Information on SNPs dataset
    paths["snps_fam"] = args.snps_matrix + ".fam"
    paths["snps_kinship"] = args.snps_matrix + ".kinship" 

    ## Creating kinship matrix specific to used accessions
    paths["pheno_intersected_fn"] = "%s/%s.phenotypes" % (args.outdir, args.name)
    paths["kinship_fn"] = "%s/%s.kinship" % (args.outdir, args.name)
    cur_cmd = "python2.7 %s --pheno %s --fam_file %s --kinship_file %s --output_pheno  %s " + \
            "--output_kinship %s --DBs_list %s"
    cur_cmd = cur_cmd % (paths["kinship_intersect_script"], paths["pheno_orig_fn"], paths["snps_fam"], 
            paths["snps_kinship"], paths["pheno_intersected_fn"], paths["kinship_fn"], args.kmers_table+".names")
    run_and_log_command(cur_cmd, f_log)
    
    ## Transformation of the phenotypes and creating permutations of phenotype
    paths["pheno_permuted_fn"] = "%s/%s.phenotypes_and_permutations" % (args.outdir, args.name)
    paths["pheno_permuted_transformed_fn"] = "%s/%s.phenotypes_permuted_transformed" % (args.outdir, args.name)
    paths["EMMA_perm_log_fn"] = "%s/EMMA_perm.log" % (args.outdir)
    paths["log_R_permute"] =    "%s/phenotypes_transformation_permutation.log" % (args.outdir)
    paths["inverse_covariance_matrix"] = "%s/inverse_covariance_matrix" % (args.outdir)
    cur_cmd = "Rscript %s %s %s %d %s %s %s %s > %s"
    cur_cmd = cur_cmd % (paths["Permute_phenotype"] , paths["pheno_intersected_fn"], paths["kinship_fn"],\
           args.n_permutations, paths["pheno_permuted_fn"], paths["pheno_permuted_transformed_fn"], 
           paths["EMMA_perm_log_fn"], paths["inverse_covariance_matrix"], paths["log_R_permute"])
    run_and_log_command(cur_cmd, f_log)
    phenotypes_names = file(paths["pheno_permuted_fn"] ,"r").read().split("\n")[0].split("\t")[1:]

    ## handles for all gemma runs (to monitor how many are running)
    gemma_handles = []
    ##############################################################################################################
    # Calculate affective minor allele frequency
    n_accession = len(file(paths["pheno_intersected_fn"]).read().split("\n")[1:-1])
    if n_accession < args.min_data_points:
        err_msg =  "Can't run with less than %d data points, there is only %d values here" % (args.min_data_points, n_accession)
        f_log.write(err_msg + "\n")
        f_log.close()
        os.system("touch %s/NOT_ENOUGH_DATA" % args.outdir)
        exit()
    affective_maf = args.maf
    maf_from_mac = float(args.mac) / float(n_accession)
    if maf_from_mac > affective_maf:
        affective_maf = maf_from_mac
    print "affective MAF =", affective_maf

    ##############################################################################################################
    ########################################      k-mers associations      #######################################
    ##############################################################################################################
    if args.run_kmers:
        paths["kmers_associations_dir"] = "%s/kmers" % args.outdir # Create relevant directory
        run_and_log_command("mkdir %s" % paths["kmers_associations_dir"], f_log)
    
        cur_cmd = "%s -p %s -b %s -o %s -n %d --parallel %d --kmers_table %s --kmer_len %d --maf %f --mac %d"  %\
                (paths["assoicte_kmers"], paths["pheno_permuted_transformed_fn"], args.name, 
                        paths["kmers_associations_dir"], args.n_kmers, args.parallel, 
                        args.kmers_table, args.kmers_len, args.maf, args.mac)
        # Optional parameters for k-mers associations
        if args.kmers_pattern_counter:
            cur_cmd = cur_cmd + " --pattern_counter"
        if args.qq_plot:
            cur_cmd = cur_cmd + " --inv_covariance_matrix " + paths["inverse_covariance_matrix"]
    
        paths["log_kmers_associations"] = "%s/associate_kmers.log" % (args.outdir)
        cur_cmd = cur_cmd + " 2> %s" % paths["log_kmers_associations"]
        run_and_log_command(cur_cmd, f_log)
    
        ## Run GEMMA on the results from the k-mers associations to get the exact scoring
        f_log.writelines("We have %d phenotypes\n" % len(phenotypes_names) )
        for (p_ind, p_name) in enumerate(phenotypes_names):
            fam_fn = get_file_type_in_dir(paths["kmers_associations_dir"], "%s.fam" % p_name)
            ## for all the fam files, move to fam.orig and create a new fam file with the non-transformed phenotype
            run_and_log_command("mv %s %s.orig" % (fam_fn, fam_fn), f_log) # Saving the original fam (which has transformed values)
            base_name = fam_fn[:-4]
            # building the fam file to run gemma
            cur_cmd = r"""cat %s | tail -n +2 | awk '{print $1 " " $1 " 0 0 0 " $%d}' > %s""" % \
                    (paths["pheno_permuted_fn"], p_ind +2,  fam_fn)
            run_and_log_command(cur_cmd, f_log)
            # run gemma on best k-mers
            cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -maf %f -miss 0.5" % \
                    (args.gemma_path, base_name,  paths["kinship_fn"], \
                    paths["kmers_associations_dir"] + "/output", p_name, affective_maf)
            run_gemma_cmd(cur_cmd, args.parallel, gemma_handles, f_log)
        
    ##############################################################################################################
    ########################################      SNPs  associations      ########################################
    ##############################################################################################################
    if args.run_one_step_snps and args.run_two_steps_snps:
        print "The program can not run both snps approximation and exact model for SNPs"
        print "You can choose up to one of this options"
        exit()

    # The two steps option is similiar to the method for k-mers, we first run GRAMMAR-Gamma for all SNPS, taking 
    # the best X and then run the full model only on this part.
    # Notice that in any case, the two steps is only for the permutation, for the real values, we run the exact
    # model on all SNPs
    if args.run_one_step_snps or args.run_two_steps_snps:
        f_log.write("\n\nStart associating snps\n")
        # Create the directory for SNPs associations
        paths["snps_associations_dir"] = "%s/snps" % args.outdir
        run_and_log_command("mkdir %s" % paths["snps_associations_dir"], f_log)
        
        # Make a link to the bed/bim of the full dataset
        paths["snps_table_fn"] = "%s/snps.plink" % paths["snps_associations_dir"]
        run_and_log_command("ln -s %s.bed %s.bed" % (args.snps_matrix, paths["snps_table_fn"]), f_log)
        run_and_log_command("ln -s %s.bim %s.bim" % (args.snps_matrix, paths["snps_table_fn"]), f_log)
        
        # building a fam file to run vs. all snps (so it will be expanded to the full format with missing values)
        create_full_fam_file("%s.fam" % paths["snps_table_fn"],paths["snps_fam"],paths["pheno_permuted_fn"],args.n_permutations+1)
        
        # Run GEMMA on the full phenotypes
        cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -n %d -maf %f -miss 0.5" % \
                (args.gemma_path, paths["snps_table_fn"], paths["snps_kinship"], \
                paths["snps_associations_dir"] + "/output", \
                phenotypes_names[0], 1, affective_maf)
        run_gemma_cmd(cur_cmd, args.parallel, gemma_handles, f_log)

        if args.run_two_steps_snps:
            cur_cmd = "%s %s %s %s %s %f %d" % \
                    (paths["associate_snps"],paths["pheno_permuted_transformed_fn"], args.snps_matrix, \
                    paths["snps_associations_dir"] + "/" + args.name, args.n_snps, affective_maf, args.mac)
            run_and_log_command(cur_cmd, f_log)

            for (p_ind, p_name) in enumerate(phenotypes_names[1:]):
                cur_base_bedbim = paths["snps_associations_dir"] + "/" + args.name + "." + p_name
                # Create the relevan fam file (just link the general one we created)
                run_and_log_command("cp %s.fam %s.fam" % (paths["snps_table_fn"], cur_base_bedbim), f_log)
                # Run GEMMA
                cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -n %d -maf %f -miss 0.5"  % \
                        (args.gemma_path, cur_base_bedbim, paths["snps_kinship"], \
                        paths["snps_associations_dir"] + "/output", \
                        p_name, p_ind+2, affective_maf)
                run_gemma_cmd(cur_cmd, args.parallel, gemma_handles, f_log)

        if args.run_one_step_snps:
            for (p_ind, p_name) in enumerate(phenotypes_names[1:]):
                cur_cmd = "%s -bfile %s -lmm 2 -k %s -outdir %s -o %s -n %d -maf %f -miss 0.5" % \
                        (args.gemma_path, paths["snps_table_fn"], paths["snps_kinship"], \
                        paths["snps_associations_dir"] + "/output", \
                        p_name, p_ind+1, affective_maf)
                run_gemma_cmd(cur_cmd, args.parallel, gemma_handles, f_log)
    
    ##############################################################################################################
    ##################################      Accumulate general statistics      ###################################
    ##############################################################################################################
    while count_running_gemma(gemma_handles) != 0: #make sure all GEMMA finished running
        time.sleep(30)
    
    ## Save the smallest p-values for k-mers associations in the permutations
    if args.run_one_step_snps or args.run_two_steps_snps:
        paths["best snps pvals"] = paths["snps_associations_dir"] + "/" + "best_pvals"
        res = calc_best_pvals(paths["snps_associations_dir"] + "/output", paths["best snps pvals"])
        th_5per = get_threshold_from_perm(res, "P", args.n_permutations, 0.05)
        th_10per = get_threshold_from_perm(res, "P", args.n_permutations, 0.1)
        run_and_log_command("echo %f > %s/threshold_5per" % (th_5per, paths["snps_associations_dir"]),f_log)
        run_and_log_command("echo %f > %s/threshold_10per" % (th_10per, paths["snps_associations_dir"]),f_log)
        run_and_log_command("cat %s/output/phenotype_value.assoc.txt | awk '(-log($9)/log(10)) > %f' > %s/pass_threshold_5per" %\
                (paths["snps_associations_dir"], th_5per, paths["snps_associations_dir"]), f_log)
        run_and_log_command("cat %s/output/phenotype_value.assoc.txt | awk '(-log($9)/log(10)) > %f' > %s/pass_threshold_10per" %\
                (paths["snps_associations_dir"], th_10per, paths["snps_associations_dir"]), f_log)
        
    if args.run_kmers:
        paths["best kmers pvals"] = paths["kmers_associations_dir"] + "/" + "best_pvals"
        res = calc_best_pvals(paths["kmers_associations_dir"] + "/output", paths["best kmers pvals"])
        th_5per = get_threshold_from_perm(res, "P", args.n_permutations, 0.05)
        th_10per = get_threshold_from_perm(res, "P", args.n_permutations, 0.1)
        run_and_log_command("echo %f > %s/threshold_5per" % (th_5per, paths["kmers_associations_dir"]),f_log)
        run_and_log_command("echo %f > %s/threshold_10per" % (th_10per, paths["kmers_associations_dir"]),f_log)
        run_and_log_command("cat %s/output/phenotype_value.assoc.txt | awk '(-log($9)/log(10)) > %f' > %s/pass_threshold_5per" %\
                (paths["kmers_associations_dir"], th_5per, paths["kmers_associations_dir"]), f_log)
        run_and_log_command("cat %s/output/phenotype_value.assoc.txt | awk '(-log($9)/log(10)) > %f' > %s/pass_threshold_10per" %\
                (paths["kmers_associations_dir"], th_10per, paths["kmers_associations_dir"]), f_log)
    
    ##############################################################################################################
    #######################################   Clean intermediate files   #########################################
    ##############################################################################################################
    if args.remove_intermediate:
        for file_type in ["bed","bim","fam"]:
            if args.run_kmers:
                run_and_log_command("rm %s/%s.*.P*.%s" % (paths["kmers_associations_dir"],args.name, file_type), f_log)
            if args.run_one_step_snps or args.run_two_steps_snps:
                run_and_log_command("rm %s/%s.P*.%s" % (paths["snps_associations_dir"],args.name, file_type), f_log)
        if args.run_one_step_snps or args.run_two_steps_snps:
            run_and_log_command("rm %s/output/P*" % paths["snps_associations_dir"], f_log)
            run_and_log_command("gzip %s/output/phenotype_value.assoc.txt"  % paths["snps_associations_dir"], f_log)
        if args.run_kmers:
            run_and_log_command("rm %s/%s.*.P*.%s" % (paths["kmers_associations_dir"],args.name, "fam.orig"), f_log)
    
    # Close log file
    f_log.close()


if __name__ == "__main__":
    main()
