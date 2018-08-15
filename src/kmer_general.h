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
#define NULL_KEY 0xFFFFFFFFFFFFFFFF

// Struct: Hash for 64 bits
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

// Define my hash tables / sets
typedef google::dense_hash_set<uint64, Hash64> kmer_set; 
typedef google::dense_hash_map<uint64, uint64, Hash64> my_hash; 

// Struct: Holds info on k-mers KMC DB
struct KMC_db_handle {
	std::string dir_path;
	std::string name;
};

// Func: return the full path of a KMD db from its handle
inline std::string KMC_db_full_path(const KMC_db_handle &h) { return h.dir_path + "/" + h.name;}

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

// Check if a kmer is part of set (can we make it one line?)
inline bool lookup_x(const kmer_set& Set, const uint64& kmer)
{
	kmer_set::const_iterator it  = Set.find(kmer);
	return (it != Set.end()); 
}

// Func: Read the list of accessions to use
std::vector<KMC_db_handle> read_accession_db_list(char *filename);
// Func: transform a bit representation to bp representation
std::string bits2kmer31(uint64 w); 
// Func: Read a file with a list of k-mers (need to add prefix/header to files)
kmer_set load_kmer_raw_file(std::string filename, std::size_t set_initial_size = 100000); 
// Func: Read a file with a list of k-mers with scores (need to add prefix/header to files)
kmer_set load_kmer_and_score_raw_file(std::string filename, std::size_t set_initial_size = 100000);

double get_time(void);

#endif

