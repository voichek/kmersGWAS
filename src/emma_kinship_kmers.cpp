///
///      @file  emma_kinship_kmers.cpp
///     @brief  calculate a kinship matrix from the presence/absence of k-mers
///
/// Using the presence / absence pattern of the k-mers we will calculate a relatedness matrix
/// using the same calculation as the EMMA.kinship one
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  12/23/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include <iostream>
#include <string>
#include <vector>

#include "kmer_general.h"
#include "kmers_multiple_databases.h"

using namespace std;

int main(int argc, char* argv[]) {
	if(argc != 4) {
		cerr << "usage: " << argv[0] << " <kmers table> <kmers_len> <MAF>" << endl;
		return -1;
	}
	string fn_kmers_table(argv[1]);
	size_t kmer_len(atoi(argv[2]));
	double MAF(atof(argv[3]));


	MultipleKmersDataBases multiDB(
			fn_kmers_table,
			load_kmers_talbe_column_names(fn_kmers_table),	
			kmer_len);

	size_t n_acc = load_kmers_talbe_column_names(fn_kmers_table).size();
	size_t min_count = ceil(static_cast<double>(n_acc) * MAF);
	cerr << "Min count = " << min_count << endl;
	uint64_t n_snps(0);
	vector<vector<uint64_t> > K(n_acc, vector<uint64_t>(n_acc, 0));
	
	cerr << "loading..." << endl;
	while(multiDB.load_kmers(1<<20, min_count)) {
		cerr << "."; cerr.flush();
		multiDB.update_emma_kinshhip_calculation(K, n_snps);
	}
	cerr << "#" << n_snps << endl;

	vector<vector<double> > K_norm(n_acc, vector<double>(n_acc, 0));
	for(size_t i=0; i<n_acc; i++) {
		K_norm[i][i] = 1;
		for(size_t j=0; j<i; j++) {
			K_norm[i][j] = static_cast<double>(K[i][j])  / static_cast<double>(n_snps);
			K_norm[j][i] = K_norm[i][j];
		}
	}

	for(size_t i=0; i<n_acc; i++) {
		for(size_t j=0; j<n_acc; j++) {
			if(j>0)
				cout << "\t";
			cout << K_norm[i][j];
		}
		cout << "\n";
	}
	return 0;
}
