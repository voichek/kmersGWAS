///
///      @file  kmers_add_strand_information.cpp
///     @brief  Take the KMC k-mers counts with and without cononization and add to the k-mers
///     on which strands there were found
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  03/10/19
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2019, Yoav Voichek
///
///This source code is released for free distribution under the terms of the
///GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmer_general.h"

#include <bitset>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

inline tuple<uint64_t, uint64_t> is_canonized_kmer_representation_flag(const uint64_t& k, const uint32& k_len) {
	uint64_t k_rc = kmer_reverse_complement(k, k_len); 
	if(k<k_rc) 
		return make_tuple(k,    0x4000000000000000);
	else
		return make_tuple(k_rc, 0x8000000000000000);
}

int main(int argc, char* argv[]) {
	if(argc != 5) {
		cerr << "usage: " << argv[0] << 
			" <KMC DB> <KMC DB non canonic kmers> <kmers len> <output file>" << endl;
		return 1;
	}
	uint32 kmer_length = atoi(argv[3]);	

	// Go over the regular KMC database kmers and add to hashtable as keys with value 0
	KmerUint64Hash kmers(1000000); // Can get the initial hash_table_size from user
	kmers.set_empty_key(NULL_KEY);


	CKmerUpTo31bpAPI kmer_obj(kmer_length); uint kmer_counter;
	CKMCFile kmer_db_canon;
	kmer_db_canon.OpenForListing(argv[1]);
	while(kmer_db_canon.ReadNextKmer(kmer_obj, kmer_counter)) 
		kmers.insert(KmerUint64Hash::value_type(kmer_obj.to_uint(),0));
	cout << "Canonized kmers:\t" << kmers.size() << endl;
	
	// Go over the non-canon' KMC database save flags (0x1 / 0x2) if apeared as canon/non-canon
	CKMCFile kmer_db_non_canon;
	kmer_db_non_canon.OpenForListing(argv[2]);
	KmerUint64Hash::iterator it_hash;
	size_t counter_all(0), counter_found(0);
	while(kmer_db_non_canon.ReadNextKmer(kmer_obj, kmer_counter)) {
		uint64_t k, flag;
		tie(k, flag)  = is_canonized_kmer_representation_flag(kmer_obj.to_uint(), kmer_length);
		it_hash = kmers.find(k);
		if(it_hash != kmers.end()) {
			it_hash->second |= flag;
			counter_found++;
		}
		counter_all++;
	}

	cout << "Non-canon kmers:\t" << counter_all << endl;
	cout << "Non-canon kmers found:\t" << counter_found << endl;
	// output to file all kmers with the flags in the most significant bits (as k-mers is of len < 31
	// there is at least two free bit
	vector<size_t> flag_counter(4,0);
	vector<uint64_t> kmers_to_save;
	for(auto it : kmers) {
		flag_counter.at(it.second >> 62)++;
		kmers_to_save.push_back(it.first | it.second); // add flags in 2-MSBs
	}
	for(size_t i=0; i<4; i++)
		cout << "flag\t" << i << "\tcount is\t" << flag_counter[i] << endl;

	sort(kmers_to_save.begin(), kmers_to_save.end(),
			[](const uint64_t & a, const uint64_t & b) -> bool
			{return (a&0x3FFFFFFFFFFFFFFF) < (b&0x3FFFFFFFFFFFFFFF);});
	cout << "kmers to save:\t" << kmers_to_save.size() << endl;

	ofstream fout(argv[4], ios::binary);
	for(size_t i=0; i<kmers_to_save.size(); i++)
		fout.write(reinterpret_cast<const char *>(&kmers_to_save[i]), sizeof(kmers_to_save[i]));
	fout.close();
	return 0;
}
