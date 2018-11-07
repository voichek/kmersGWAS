///
///      @file  F_create_kmer_table.cpp
///     @brief  This program will create a table of presence absence accross DBs
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  10/29/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///


#include "kmer_general.h"
#include "kmer_multipleDB.h"

using namespace std;

int main(int argc, char *argv[]) {
	/* Read user input */
	if(argc != 5) {
		cerr << "usage: " << argv[0] << 
			" <file with KMC DBs paths> <name of sorted k-mer file> <output file> <kmer len>" << endl;
		return -1;
	}

	// Read accessions to use
	vector<KMC_db_handle> db_handles = read_accession_db_list(argv[1]); 
	vector<string> db_paths, db_names;
	
	for(size_t i=0; i<db_handles.size(); i++) {
		db_names.push_back(db_handles[i].name);
		db_paths.push_back(db_handles[i].dir_path);
	}
	cerr << "Create merger" << endl;	
	kmer_multipleDB_merger merger(db_paths, db_names, argv[2], atoi(argv[4]));
	
	cerr << "open file" << endl;
	ofstream fout(argv[3], ios::binary);
	merger.output_table_header(fout);	
	uint64_t total_iter = 5000;
	for(uint i=1; i<=total_iter; i++) {
		cerr << i << " - " << total_iter << " : Loading k-mers" << endl;
		merger.load_kmers(i, total_iter);
		cerr << "Outputing to table" << endl;
		merger.output_to_table(fout);
	}	
	cerr << "close file" << endl;
	fout.close();
	return 0;
}
