///
///      @file  histogram_KMC_kmers_counts.cpp 
///     @brief  Create an histogram of k-mers count for KMC DB 
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

#include "kmers_single_database.h"

using namespace std;
int main(int argc, char *argv[]) {
	if(argc != 4) {
		cerr << "usage: " << argv[0] << " <path of DB> <DB name> <kmer length>" << endl;
		return 1;
	}
	string dir_path(argv[1]);
	string db_name(argv[2]);
	KmersSingleDataBase acc_kmers(dir_path, db_name, atoi(argv[3]));

	vector<size_t> counters = acc_kmers.calculate_kmers_counts_histogram();
	
	cout << "appearance\tcount" << endl;
	for(size_t i=0; i<counters.size(); i++) 
		cout << i << "\t" << counters[i] << endl;											
	
	return 0;
}

