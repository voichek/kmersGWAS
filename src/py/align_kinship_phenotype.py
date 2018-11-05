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
import sys
import argparse
## General parameters:

execfile("%s/functions.py" % sys.path[0])


parser = argparse.ArgumentParser(description="Create a sub-matrix only of samples found in the phenotype file")
parser.add_argument("--pheno", dest = "fn_phenotype", type=str, required=True,
        help='phenotype file name Format:sample name[TAB]phenotype val[NEW-LINE]...)')

parser.add_argument("--fam_file", dest = "fn_fam", type=str, required=True,
        help="Fam file with the order of the samples")

parser.add_argument("--kinship_file", dest = "fn_kinship", type=str, required=True,
        help="file with the kinship matrix")

parser.add_argument("--output_pheno", dest = "fn_out_pheno", type=str, default="",
        help="Name of output phenotype file")

parser.add_argument("--output_kinship", dest = "fn_out_kinship", type=str, required=True,
        help="Name of output kinship file")

parser.add_argument("--DBs_list", dest = "db_list", type=str, required=False, default = "",
        help="list of DBs, if provided output phenotypes/kinship matrix will be intersected with it")

parser.add_argument("-v", "--verbose", help="increase output verbosity",
                action="store_true")



def main():
    args = parser.parse_args()

    # Create a list of accessions found in the kinship matrix and in the phenotype file
    kinship_acc = get_column(args.fn_fam, 0, " ")
    phenotype_acc = get_column(args.fn_phenotype, 0)[1:]
    phenotype_val = get_column(args.fn_phenotype, 1)[1:] # Phenotypes values
    intersect_acc = [x for x in phenotype_acc if x in kinship_acc]
    if args.verbose:
        print "We have %d accession in the kinship matrix" % len(kinship_acc)
        print "We have %d in this phenotype" % len(phenotype_acc)
        print "We have %d in the intersection" % len(intersect_acc)
    
    if len(args.db_list) > 0: 
        intersect_acc = [x for x in intersect_acc if x in get_column(args.db_list, 1)]
        if args.verbose:
            print "Intersection also with external list (DBs), left with %d" % len(intersect_acc)

    index_kinship = [] ## Finding the indices of kinship in the new list 
    for i in range(len(intersect_acc)):
        index_kinship.append(kinship_acc.index(intersect_acc[i]))
    create_new_relatedness_matrix(args.fn_kinship, args.fn_out_kinship, index_kinship)
    
    if len(args.fn_out_pheno) > 0:
        # Create new phenotype file
        fout_pheno = file(args.fn_out_pheno, "w")
        fout_pheno.write("accession_id\tphenotype_value\n")
        for i in range(len(intersect_acc)):
            cur_ind = phenotype_acc.index(intersect_acc[i])
            fout_pheno.write("%s\t%s\n" % (phenotype_acc[cur_ind], phenotype_val[cur_ind]))
        fout_pheno.close()
    
if __name__ == "__main__":
    main()
