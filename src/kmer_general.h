#ifndef KMER_GENERAL_H
#define KMER_GENERAL_H

#include <sys/time.h>
#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <sparsehash/dense_hash_map>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/dense_hash_set>

#include "KMC/kmc_api/kmc_file.h"

#define MAX_KMER_LEN 31 
#define MIN_KMER_LEN 15 //Arbitrary

#define WLEN 64
#define NULL_KEY 0xFFFFFFFFFFFFFFFF

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

// Definition of hash tables / sets
typedef google::dense_hash_set<uint64_t, Hash64> KmersSet; 
typedef google::dense_hash_map<uint64_t, uint64_t, Hash64> KmerUint64Hash; 
typedef google::sparse_hash_map<uint64_t, uint64_t, Hash64> KmerUint64SparseHash; 

// Struct: Holds info on k-mers KMC DB
struct AccessionPath {
	std::string path;
	std::string name;
};

typedef std::pair <std::vector<std::string>, std::vector<float> > PhenotypeList;

std::vector<std::string> load_kmers_talbe_column_names(const std::string &kmers_table_base);

// Func: return the full path of a KMD db from its handle
inline std::string KMC_db_full_path(const AccessionPath &h) { return h.path + "/" + h.name;}

void filter_kmers_to_set(std::vector<uint64_t> &kmers, const KmersSet &set_kmers);

class CKmerUpTo31bpAPI: public CKmerAPI {
	public:
		CKmerUpTo31bpAPI (uint32 length = 0): CKmerAPI(length), 
		m_shift(64 - (((kmer_length - 1 + byte_alignment) % 32) * 2) -2) {
			if(length>31)
				throw std::invalid_argument("k-mer length should be <=31");
		}
		uint64_t to_uint() {return (uint64_t)kmer_data[0] >> m_shift;}
		uint32 get_kmer_length() const {return kmer_length;}
		void plot_info() {
			std::cerr << "byte_alignment:\t" << (uint64_t)byte_alignment << std::endl;
			std::cerr << "no_of_rows:\t"    << no_of_rows << std::endl;
			std::cerr << "shift:\t" << m_shift << std::endl;
		}
	private:
		uint32 m_shift;		
};

// Check if a kmer is part of set (can we make it one line?)
inline bool lookup_x(const KmersSet& Set, const uint64_t& kmer)
{
	KmersSet::const_iterator it  = Set.find(kmer);
	return (it != Set.end()); 
}

// Func: Read the list of accessions to use
std::vector<AccessionPath> read_accessions_path_list(std::string filename);
// Func: transform a bit representation to bp representation
std::string bits2kmer31(uint64_t w, const std::size_t& k); 
// Func: Read a file with a list of k-mers (need to add prefix/header to files)
KmersSet load_kmer_raw_file(std::string filename, 
		std::size_t set_initial_size = 100000, 
		const bool with_scores = false); 
// Func: Read a file with a list of k-mers with scores (need to add prefix/header to files)
KmersSet load_kmer_and_score_raw_file(std::string filename, std::size_t set_initial_size = 100000);
// Func: bitwise-reverse of uint64
uint64_t reverseOne(uint64_t x);

// Func: kmer reverse complement
inline uint64_t kmer_reverse_complement(uint64_t x, const uint32& k_len) {
	x = ((x & 0xFFFFFFFF00000000) >> 32) | ((x & 0x00000000FFFFFFFF) << 32);
	x = ((x & 0xFFFF0000FFFF0000) >> 16) | ((x & 0x0000FFFF0000FFFF) << 16);
	x = ((x & 0xFF00FF00FF00FF00) >> 8)  | ((x & 0x00FF00FF00FF00FF) << 8);
	x = ((x & 0xF0F0F0F0F0F0F0F0) >> 4)  | ((x & 0x0F0F0F0F0F0F0F0F) << 4);
	x = ((x & 0xCCCCCCCCCCCCCCCC) >> 2)  | ((x & 0x3333333333333333) << 2);
	return (~x) >> (64 - k_len - k_len);
}

double get_time(void);

typedef std::tuple<uint64_t, double, size_t> AssociationScoreHeap; //k-mer / score/ row index in table
typedef std::tuple<uint64_t, uint64_t, size_t> AssociationOutputInfo; // k-mer name / row index in table

struct kmers_output_list
{
	std::vector<AssociationOutputInfo> list;
	size_t next_index;
};

struct cmp_second
{
	inline bool operator() (const AssociationScoreHeap& left, const AssociationScoreHeap& right) const
	{return (std::get<1>(left)) > (std::get<1>(right));} 
};

typedef std::priority_queue<AssociationScoreHeap, std::vector<AssociationScoreHeap>, cmp_second> AssociationsPriorityQueue;

/**
 * @struct BedBimFilesHandle
 * @brief  holds handles to a bed & bim files (and write the header of the bed file)
 */
struct BedBimFilesHandle {
	BedBimFilesHandle(const std::string &base_name):
		f_bed(base_name + ".bed",std::ios::binary),
		f_bim(base_name + ".bim", std::ios::out) 
	{ f_bed << (char)0x6C << (char)(0x1B) << (char)(0x01);} // header of bed file
	~BedBimFilesHandle() {close();}
	BedBimFilesHandle(BedBimFilesHandle&& o) = default;
	void close();

	std::ofstream f_bed;
	std::ofstream f_bim;
};

// Permuting the order of the phenotypes for the implementation of the scoring
void permute_scores(std::vector<float> &V);  // assume V is a multiplication of 128

std::pair<std::vector<std::string>, std::vector<PhenotypeList>> load_phenotypes_file(const std::string &filename);

void write_fam_file(const std::vector<PhenotypeList> &phenotypes, const std::string &fn);

void write_fam_file(const PhenotypeList &phenotype, const std::string &fn);

std::size_t get_index_DB(const std::string &name, const std::vector<std::string> &names);

PhenotypeList intersect_phenotypes_to_present_DBs(const PhenotypeList &pl, const std::string &kmers_table_base, const bool &must_be_present);

uint64_t kmers_step_to_threshold(const uint64_t &step, const uint64_t &total_steps, const uint64_t &kmer_length);

uint64_t kmer2bits(std::string k);

bool is_file_exist(const char *fileName);
bool is_file_exist(const std::string &filneName);

#endif

