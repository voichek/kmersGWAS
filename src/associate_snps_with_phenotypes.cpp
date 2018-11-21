///
///      @file  associate_snps_with_phenotypes.cpp
///     @brief  Associate snps with phenotypes and output the best n snps to check with a 
//	non approximated method
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  11/20/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///


#include "kmer_general.h"
#include "snps_multiple_databases.h"
#include "best_associations_heap.h"

#include <sys/sysinfo.h> // To monitor memory usage
#include <algorithm>    // 
#include <iostream>
#include <iterator>
#include <utility> //std::pair


using namespace std;

int main(int argc, char* argv[])
{
	if(argc != 6) {
		cerr << "usage: " << argv[0] << " <phenotypes file> <base bedbim file> <base output files> <snps to output> <maf>" << endl;
		return 1;
	}
	// load phenotypes
		
	pair<vector<string>, vector<PhenotypeList> > phenotypes_info = load_phenotypes_file(argv[1]);
	cerr << "Loading snps information" << endl;
	MultipleSNPsDataBases snps_dataset(argv[2], phenotypes_info.second[0].first); // base bedbim, names of the phenotyped accessions

	size_t n_samples = phenotypes_info.second[0].first.size();
	size_t n_best_snps_to_save = atoi(argv[4]);
	
	double maf = atof(argv[5]); //minor allele frequency
	cerr << "MAF = " << maf << " n_sample = " << n_samples << endl;
	double mac = ceil(maf*n_samples); //minor allele count
	cerr << "MAC = " << mac << endl;
	size_t phenotype_n = phenotypes_info.first.size();
	vector<vector<size_t> > best_snps_indices;
	cerr << "Associating phenotypes:"; cerr.flush();
	double t0 = get_time();
	for(size_t phenotype_i=0; phenotype_i<phenotype_n; phenotype_i++) {
		cerr << "."; cerr.flush();
		best_snps_indices.push_back(snps_dataset.get_most_associated_snps(
					phenotypes_info.second[phenotype_i].second,n_best_snps_to_save, mac));

	}
	double d_time = (get_time() - t0) / phenotype_n;
	cerr << "Average time per phenotype:\t"<< d_time << endl;
	cerr << "\noutputting best snps";
	vector<string> output_base_filenames;
	for(size_t phenotype_i=0; phenotype_i<phenotype_n; phenotype_i++) 
		output_base_filenames.push_back(string(argv[3]) + "." + phenotypes_info.first[phenotype_i]);
	snps_dataset.output_plink_bed_file(output_base_filenames, best_snps_indices);
	return 0;
}
