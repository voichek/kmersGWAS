/**
 *       @file  F_kmers_intersect_and_sort.cpp
 *      @brief  Loading a list of k-mers and intersecting from each given KMC DB only these 
 *      k-mers (and plotting them in a sorted order)
 *
 * Detailed description starts here.
 *
 *     @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
 *
 *   @internal
 *     Created  08/14/18
 *    Compiler  gcc/g++
 *     Company  Max Planck Institute for Developmental Biology Dep 6
 *   Copyright  Copyright (c) 2018, Yoav Voichek
 *
 * This source code is released for free distribution under the terms of the
 * GNU General Public License as published by the Free Software Foundation.
 * =====================================================================================
 */

#include "kmer_general.h"
#include "kmers_single_database.h"

using namespace std;

int main(int argc, char *argv[]) {
	/* Read user input */
	if(argc != 5) {
		cerr << "usage: " << argv[0] << 
			" <file with KMC DBs paths> <kmers to intersect file> <output base name> <kmer len>" << endl;
		return -1;
	}

	cerr << "File with DBs filenames: " << argv[1] << endl;
	cerr << "Kmers to intersect with file: " << argv[2] << endl;
	string base_output_filename(argv[3]);
	cerr << "Output filename: " << base_output_filename << endl;

	/* Load user input files */	
	vector<KMC_db_handle> db_handles = read_accession_db_list(argv[1]); // Read accessions to use
	cerr << "Loading k-mers file" << endl;
	kmer_set kmer_list_to_use =  load_kmer_raw_file(argv[2]); // loading kmers file
	cerr << "Hash set size = " << kmer_list_to_use.size() << endl;
	
	/* Go over all input DBs and intersect k-mers with given list */
	for(uint64_t i=0; i<db_handles.size(); i++) {
		cerr << i << ". Opening DB: " << db_handles[i].name << endl;
		kmer_DB cur_acc(db_handles[i].dir_path, db_handles[i].name, atoi(argv[4]));
		cur_acc.intersect_kmers(kmer_list_to_use, base_output_filename);
	}
	return 0;
}
