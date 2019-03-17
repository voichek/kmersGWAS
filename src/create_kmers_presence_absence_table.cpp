///
///      @file  create_kmers_presence_absence_table.cpp 
///     @brief  This program will create a table of presence absence accross DBs
//		
//		Taking the sorted general kmers list as well as the sorted kmers from 
//		each accessions, create a table with kmers as rows and accessions as collumns
//		where each cell indicate if the kmers was found in the accessions or not
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
#include "kmers_merge_multiple_databaes.h"

using namespace std;

int main(int argc, char *argv[]) {
	/* Read user input */
	if(argc != 5) {
		cerr << "usage: " << argv[0] << 
			" <1 list of kmers files> <2 all sorted kmers> <3 output name> <4 kmer len>" << endl;
		return -1;
	}

	// Read accessions to use
	vector<AccessionPath> kmers_handles = read_accessions_path_list(argv[1]); 
	vector<string> kmers_filenames, accessions_names;
	
	ofstream fout_names(string(argv[3]) + ".names", ios::binary);
	for(size_t i=0; i<kmers_handles.size(); i++) {
		accessions_names.push_back(kmers_handles[i].name);
		fout_names << kmers_handles[i].name << endl;
		kmers_filenames.push_back(kmers_handles[i].path);
	}
	fout_names.close();
	cerr << "Create merger" << endl;	
	MultipleKmersDataBasesMerger merger(kmers_filenames, accessions_names, argv[2], atoi(argv[4]));
	
	cerr << "Opens file" << endl;
	ofstream fout(string(argv[3]) + ".table", ios::binary);
	merger.output_table_header(fout);	
	uint64_t total_iter = 5000;
	for(uint i=1; i<=(total_iter+1); i++) { // +1 is for debugging - should be empty
		cerr << i << " / " << total_iter << " : Loading k-mers" << endl;
		merger.load_kmers(i, total_iter);
		merger.output_to_table(fout);
	}	
	cerr << "close file" << endl;
	fout.close();
	return 0;
}
