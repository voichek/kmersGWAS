#ifndef KMER_GENERAL_H
#define KMER_GENERAL_H

#include <sys/time.h>

#include <iostream>
#include <string>
#include <vector>
#include <bitset>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

#include "../include/KMC/kmc_api/kmc_file.h"

#define KMER_LEN 31
#define WLEN 64
#define HASH_TABLE_SIZE 1300000000


struct Hash64 {
	size_t operator()(uint64_t key) const { 
		key ^= key >> 33;
		key *= 0xff51afd7ed558ccd;
		key ^= key >> 33;
		key *= 0xc4ceb9fe1a85ec53;
		key ^= key >> 33;
		return key;
	}
};

typedef google::dense_hash_set<uint64, Hash64> kmer_set; 
typedef google::dense_hash_map<uint64, uint64, Hash64> my_hash; 

// Func: Read the list of accessions to use
inline std::vector<std::string> read_accession_db_list(char *filename) {
	std::ifstream fin(filename);
	std::vector<std::string> res;
	std::string name;
	while(fin >> name) {res.push_back(name);}
	return res;
}

///
/// @brief  transform a bit representation to bp representation
/// @param  
/// @return 
///
inline std::string bits2kmer31(uint64 w) {
	const static char dict_bp[] = {'A','C','G','T'};
	const static uint64 mask2bits = 0x0000000000000003;

	std::string res(31,'X');
	for(std::size_t i=0; i<31; i++) {
		res[30-i] = dict_bp[w & mask2bits];
		w = (w>>2);
	}
	return res;
}

inline double get_time(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + (tv.tv_usec / 1000000.0);
}

class CKmerAPI_YV: public CKmerAPI {
	public:
		CKmerAPI_YV (uint32 length = 0): CKmerAPI(length) {}
		uint64 to_uint() {return (uint64)kmer_data[0];}
		void infoYV() {
			cerr << "kmer_length = " <<  kmer_length << endl;				// Kmer's length, in symbols
			cerr << "byte_alignment = " <<  (int)byte_alignment << endl;			// A number of "empty" symbols placed before prefix to let sufix's symbols to start with a border of a byte
			cerr << "no_of_rows = " <<  no_of_rows << endl;				// A number of 64-bits words allocated for kmer_data 	
		}
};


inline bool lookup_x(const kmer_set& Set, const uint64& kmer)
{
	kmer_set::const_iterator it  = Set.find(kmer);
	return (it != Set.end()); 
}

inline kmer_set load_kmer_raw_file(std::string filename, std::size_t set_initial_size = 100000) {
	kmer_set kmer_list_to_use(set_initial_size);
	kmer_list_to_use.set_empty_key(-1); // need to define empty value for google dense hash table

	ifstream kmer_file(filename, std::ifstream::binary);
	if(kmer_file) { // if file could be open
		kmer_file.seekg(0, kmer_file.end); // 
		uint64 kmers_in_file = (kmer_file.tellg()) >> 3;
		kmer_file.seekg(0, kmer_file.beg);
		uint64 kmer_uint;
		for(uint64 i=0; i<(kmers_in_file); i++) {
			kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
			kmer_list_to_use.insert(kmer_uint);
		}
		kmer_file.close();
	}
	return kmer_list_to_use;
}

inline kmer_set load_kmer_and_score_raw_file(std::string filename, std::size_t set_initial_size = 100000) {
	kmer_set kmer_list_to_use(set_initial_size);
	kmer_list_to_use.set_empty_key(-1); // need to define empty value for google dense hash table

	ifstream kmer_file(filename, std::ifstream::binary);
	if(kmer_file) { // if file could be open
		kmer_file.seekg(0, kmer_file.end); // 
		uint64 kmers_in_file = (kmer_file.tellg()) >> (3+1);
		kmer_file.seekg(0, kmer_file.beg);
		uint64 kmer_uint;
		double score;
		for(uint64 i=0; i<(kmers_in_file); i++) {
			kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
			kmer_file.read(reinterpret_cast<char *>(&score), sizeof(score));
			kmer_list_to_use.insert(kmer_uint);
		}
		kmer_file.close();
	}
	std::cerr << "loaded set of k-mers #" << kmer_list_to_use.size() << std::endl;
	return kmer_list_to_use;
}

#endif

