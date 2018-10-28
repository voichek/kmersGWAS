///
///      @file  F_kmers_count_histogram.cpp
///     @brief  Create an histogram of k-mers count for KMC DB 
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  08/14/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmer_DB.h"
int main(int argc, char *argv[]) {
	if(argc != 3) {
		cerr << "usage: " << argv[0] << " <path of DB> <DB name>" << endl;
		return 1;
	}
	string dir_path(argv[1]);
	string db_name(argv[2]);
	kmer_DB acc_kmers(dir_path, db_name);

	std::vector<std::size_t> counters = acc_kmers.calculate_kmers_counts_histogram();
	
	cout << "appearance\tcount" << endl;
	for(std::size_t i=0; i<counters.size(); i++) {
		cout << i << "\t" << counters[i] << endl;											
	}

	return 0;
}
//}

