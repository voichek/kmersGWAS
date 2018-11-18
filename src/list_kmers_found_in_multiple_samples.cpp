///
///      @file  F_list_kmers_found_in_multiple_DBs.cpp
///     @brief  Go over a list of k-mer DBs and output only k-mers appearing in at least 
/// 		N (defined by user) DBs. Then outputing the k-mers in a binary format.
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  08/14/18
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

using namespace std;

int main(int argc, char *argv[]) {
	/* Read user input */
	if(argc != 6) {
		cerr << "usage: " << argv[0] << 
			" <1 file with KMC DBs paths> <2 output file> <3 minimum k-mer counts> <4 kmer len> <5  hash table initial size>" << endl;
		return -1;
	}
	// Read accessions to use
	vector<KMC_db_handle> db_handles = read_accession_db_list(argv[1]); 
	string output_fn(argv[2]);
	size_t minimum_kmer_count = atoi(argv[3]);

	/* Build the hash table that will contain all the kmers counts */
	my_hash main_db(atoi(argv[5])); // Can get the initial hash_table_size from user
	main_db.set_empty_key(NULL_KEY); // need to define "empty key" 
	
	/* Defining variables to use while going over the k-mers */
	cerr << "k-mer length " << atoi(argv[4]) << endl;
	CKmerAPI_YV kmer_obj(atoi(argv[4]));
	vector<uint64_t> k_mers;
	uint kmer_counter;
	my_hash::iterator it_hash;

	/* Going over all k-mers DBs */
	for(size_t i=0; i<db_handles.size(); i++) {
		cerr << i << "\tloading: " << db_handles[i].name << endl;

		CKMCFile kmer_db;
		kmer_db.OpenForListing(KMC_db_full_path(db_handles[i])); // Open a KMC DB
		k_mers.resize(0);  	
		
		while (kmer_db.ReadNextKmer(kmer_obj, kmer_counter))// Reading k-mers in file
			k_mers.push_back(kmer_obj.to_uint());

		/* Update hash table */
		uint64_t or_all = 0;
		for(size_t k=0; k<k_mers.size(); k++) {
			it_hash = main_db.find(k_mers[k]);
			or_all |= k_mers[k];
			if(it_hash == main_db.end()){ 
				main_db.insert(my_hash::value_type(k_mers[k], 1));
			} else { 
				it_hash->second++;	}
		}
		cerr << bitset<64>(or_all) << endl;
	}
	/* sort all k-mers */
	k_mers.resize(0);
	vector<uint64_t> shareness(db_handles.size()+1, 0);
	for(auto it : main_db) {
		if(it.second > db_handles.size()) {
			cout << "Error! " << bitset<64>(it.first) << "\t" << it.second << endl;
		} else {
		shareness[it.second]++;
		}
		if(it.second>=minimum_kmer_count) 
			k_mers.push_back(it.first);	
	}
	cerr << "Sorting..." << endl;
	sort(k_mers.begin(), k_mers.end());	
	cerr << "Finish sorting" << endl;

	/* output all the k-mers apearing more than once to a file */
	ofstream file(output_fn, ios::binary);
	for(size_t i=0; i<k_mers.size(); i++) 
		file.write(reinterpret_cast<const char *>(&k_mers[i]), sizeof(k_mers[i]));
	file.close();

	/* Output shareness measure */
	cout << "k-mer appearance\tcount" << endl;
	for(size_t i=0; i<shareness.size(); i++) 
		cout << i << "\t" << shareness[i] << endl;

	return 0;
}
