## 
##       @file  align_kinship_phenotype.py
##      @brief  Given phenotype list and a kinship matrix, we will take the intersect between 
##              the accessions found in both, and order them in the same way
## Detailed description starts here.
## 
##     @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
## 
##   @internal
##     Created  09/23/18
##    Revision  $Id: doxygen.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
##     Company  Max Planck Institute for Developmental Biology Dep 6
##   Copyright  Copyright (c) 2018, Yoav Voichek
## 
## This source code is released for free distribution under the terms of the
## GNU General Public License as published by the Free Software Foundation.
## =====================================================================================
## 

## General parameters:
base_path = "/ebio/abt6/yvoichek/1001G_1001T_comparison/code/k_mer_clusters/acc_kmer_counts/correlate_phenotype/"
execfile(base_path + "scripts/general_params.py")

def main():
    ## Specific run parameters:
    fn_phenotype = sys.argv[1]
    output_basename = sys.argv[2]
    if len(sys.argv) == 4: ## We also have additional list to intersect (DB list)
        additional_acc_list = get_column(sys.argv[3], 1)
    else:
        additional_acc_list = []

    # new file for kinship matrix
    fn_new_kn_gemma = output_basename + ".kinship.gemma"
    fn_new_kn_emma = output_basename + ".kinship.emma"
    # new file for phenotypes
    fn_new_pheno = output_basename + ".pheno"
    
    # Create a list of accessions found in the kinship matrix and in the phenotype file
    kinship_acc = get_column(original_fam, 0, " ")
    phenotype_acc = get_column(fn_phenotype, 0)[1:]
    phenotype_val = get_column(fn_phenotype, 1)[1:] # Phenotypes values
    intersect_acc = [x for x in phenotype_acc if x in kinship_acc]
    print "We have %d accession in the kinship matrix" % len(kinship_acc)
    print "We have %d in this phenotype" % len(phenotype_acc)
    print "We have %d in the intersection" % len(intersect_acc)
    if len(additional_acc_list) > 0:
        intersect_acc = [x for x in intersect_acc if x in additional_acc_list]
        print "Intersection also with external list (DBs), left with %d" % len(intersect_acc)
    # We want both kinship and phenotype file to be in the order of intersect_acc

    index_kinship = [] ## Finding the indices of kinship in the new list 
    for i in range(len(intersect_acc)):
        index_kinship.append(kinship_acc.index(intersect_acc[i]))
    create_new_relatedness_matrix(kinship_matrix_EMMA, fn_new_kn_emma, index_kinship)
    create_new_relatedness_matrix(kinship_matrix_GEMMA, fn_new_kn_gemma, index_kinship)
    
    # Create new phenotype file
    fout_pheno = file(fn_new_pheno, "w")
    fout_pheno.write("accession_id\tphenotype_value\n")
    for i in range(len(intersect_acc)):
        cur_ind = phenotype_acc.index(intersect_acc[i])
        fout_pheno.write("%s\t%s\n" % (phenotype_acc[cur_ind], phenotype_val[cur_ind]))
    fout_pheno.close()
    
if __name__ == "__main__":
    main()
