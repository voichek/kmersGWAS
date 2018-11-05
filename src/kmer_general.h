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

#define MAX_KMER_LEN 31 
#define MIN_KMER_LEN 15 //Arbitrary

#define WLEN 64
#define NULL_KEY 0xFFFFFFFFFFFFFFFF

// Struct: Hash for 64 bits

/**
 * @brief   Hash function for 64 bit
 * @param   get a 64 bit unsigned long long
 * @return  size_t (32 bit)
 */
struct Hash64 {
	std::size_t operator()(uint64_t key) const { 
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

typedef std::pair <std::vector<std::string>, std::vector<double> > phenotype_list;
typedef std::pair <std::vector<std::string>, std::vector<uint64> > uint64_phenotype_list;

// Func: return the full path of a KMD db from its handle
inline std::string KMC_db_full_path(const KMC_db_handle &h) { return h.dir_path + "/" + h.name;}

void filter_kmers_to_set(std::vector<uint64> &kmers, const kmer_set &set_kmers);

class CKmerAPI_YV: public CKmerAPI {
	public:
		CKmerAPI_YV (uint32 length = 0): CKmerAPI(length), 
		m_shift(64 - (((kmer_length - 1 + byte_alignment) % 32) * 2) -2) {
			if(length>31)
				throw std::invalid_argument("k-mer length should be <=31");
		}
		uint64 to_uint() {return (uint64)kmer_data[0] >> m_shift;}
		uint32 get_kmer_length() const {return kmer_length;}
		void plot_info() {
			std::cerr << "byte_alignment:\t" << (uint64)byte_alignment << std::endl;
			std::cerr << "no_of_rows:\t"    << no_of_rows << std::endl;
			std::cerr << "shift:\t" << m_shift << std::endl;
		}
	private:
		uint32 m_shift;		
};

// Check if a kmer is part of set (can we make it one line?)
inline bool lookup_x(const kmer_set& Set, const uint64& kmer)
{
	kmer_set::const_iterator it  = Set.find(kmer);
	return (it != Set.end()); 
}

// Func: Read the list of accessions to use
std::vector<KMC_db_handle> read_accession_db_list(std::string filename);
// Func: transform a bit representation to bp representation
std::string bits2kmer31(uint64 w, const std::size_t& k); 
// Func: Read a file with a list of k-mers (need to add prefix/header to files)
kmer_set load_kmer_raw_file(std::string filename, std::size_t set_initial_size = 100000, const bool with_scores = false); 
// Func: Read a file with a list of k-mers with scores (need to add prefix/header to files)
kmer_set load_kmer_and_score_raw_file(std::string filename, std::size_t set_initial_size = 100000);
// Func: bitwise-reverse of uint64
uint64 reverseOne(uint64 x);
double get_time(void);

#endif

