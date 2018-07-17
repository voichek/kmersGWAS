#ifndef KMER_GENERAL_H
#define KMER_GENERAL_H

#include <iostream>
#include <string>
#include <vector>
#include <sys/time.h>
#include <bitset>
#include <fstream>
#include <unordered_map>
#include <stdlib.h>
#include <algorithm>
#include <cstdlib>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include "../include/KMC/kmc_api/kmc_file.h"
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>

#define KMER_LEN 31
#define WLEN 64


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

// Func: Read the list of accessions to use
inline std::vector<std::string> read_accession_db_list(char *filename) {
	std::ifstream fin(filename);
	std::vector<std::string> res;
	std::string name;
	while(fin >> name) {res.push_back(name);}
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
		uint64 to_uint() {
			return (uint64)kmer_data[0];
		}
		void infoYV() {
			cerr << "kmer_length = " <<  kmer_length << endl;				// Kmer's length, in symbols
			cerr << "byte_alignment = " <<  (int)byte_alignment << endl;			// A number of "empty" symbols placed before prefix to let sufix's symbols to start with a border of a byte
			cerr << "no_of_rows = " <<  no_of_rows << endl;				// A number of 64-bits words allocated for kmer_data 	
		}
};


#endif

