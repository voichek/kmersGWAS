
#include "kmer_general.h"
#include "kmers_single_database.h"
#include <fstream>

using namespace std;

KmersSingleDataBase::KmersSingleDataBase(const string& dir_path, const string& db_name, const uint32& kmer_length):
	m_db_name(db_name),
	m_dir_path(dir_path),
	m_sorted_kmers_f(),
	m_kmer_len(kmer_length)
{
	if((m_kmer_len<MIN_KMER_LEN) || (m_kmer_len>MAX_KMER_LEN)) {
			throw std::logic_error("Ilegal k-mer length");
	}
}


// go over the KMC DB and output to file only the kmers in the kmers set inputed
void KmersSingleDataBase::intersect_kmers(const KmersSet& kmers_to_use, std::string file_name) {
	string output_file = m_dir_path + "/" +  file_name;
	ofstream of(output_file, ios::binary);

	CKMCFile kmer_database = get_KMC_handle(); //open DB
	CKmerUpTo31bpAPI kmer_obj(m_kmer_len);

	unsigned int counter;
	uint64_t cnt_f = 0, cnt_nf = 0, kmer;

	vector<uint64_t> kmers;
	while(kmer_database.ReadNextKmer(kmer_obj, counter)) {
		kmer = kmer_obj.to_uint();
		if(lookup_x(kmers_to_use, kmer)) {
			cnt_f++;
			kmers.push_back(kmer);
		} else { cnt_nf++;}
	}
	// Sort kmers 
	sort(kmers.begin(), kmers.end());
	for(vector<uint64_t>::const_iterator it = kmers.begin(); it != kmers.end(); ++it) 
		of.write(reinterpret_cast<const char *>(&(*it)), sizeof(*it));

	of.close();
}

// Counts how many times each k-mer appeared and plot to std::cout
vector<size_t> KmersSingleDataBase::calculate_kmers_counts_histogram() {
	CKMCFile kmer_database = get_KMC_handle();
	CKmerUpTo31bpAPI kmer_obj(m_kmer_len);
	vector<size_t> counters(0, 0);
	unsigned int counter;
	while(kmer_database.ReadNextKmer(kmer_obj, counter)) {
		if(counter >= counters.size()) 
			counters.resize(counter+1,0);
		counters[counter] += 1;
	}
	return counters;
}

// get a hanle to the KMC DB to run over its k-mers 
CKMCFile KmersSingleDataBase::get_KMC_handle() {
	CKMCFile kmer_database;
	kmer_database.OpenForListing(m_dir_path + "/" +  m_db_name);

	return kmer_database;
}

/**
 * Definition of class KmersSingleDataBaseSortedFile member functions
 */
KmersSingleDataBaseSortedFile::KmersSingleDataBaseSortedFile():
	m_fin(),
	m_last_kmer(NULL_KEY),
	m_flag(NULL_KEY),
	m_kmers_in_file(NULL_KEY),
	m_kmers_count(NULL_KEY) {
	}


KmersSingleDataBaseSortedFile::KmersSingleDataBaseSortedFile(const std::string &filename):
	KmersSingleDataBaseSortedFile() {
		open_file(filename);
	}

KmersSingleDataBaseSortedFile::~KmersSingleDataBaseSortedFile() {// Make sure the file is close
	if(m_fin.is_open()) {
		close_file();
	}
}
void KmersSingleDataBaseSortedFile::open_file(const std::string &filename) {
	if(m_fin.is_open()) { // First check if file is allready open
		close_file();
	}

	m_fin.open(filename, ios::binary); // opening new file
	if(m_fin.is_open()) {
		/* Counting the number of kmers in file */
		m_fin.seekg(0, m_fin.end); // 
		m_kmers_in_file = (m_fin.tellg()) >> 3; // divide by 8 as we count bytes (kmer = 8b)
		m_fin.seekg(0, m_fin.beg);

		if(m_kmers_in_file>0) {
			m_kmers_count = 0; // First k-mer (index from 1..)
			read_kmer();
		}
		else {
			throw std::logic_error("sorted kmer file is empty: " + filename);
		}
	} 
	else {
		throw std::logic_error("can't open file: " + filename);		
	}
}

void KmersSingleDataBaseSortedFile::close_file() {
	if(m_fin.is_open()) {
		m_fin.close();
		m_flag = NULL_KEY;
		m_last_kmer = NULL_KEY;
		m_kmers_in_file = NULL_KEY;
		m_kmers_count = NULL_KEY;
	}
}

void KmersSingleDataBaseSortedFile::read_kmer() {
	m_fin.read(reinterpret_cast<char *>(&m_last_kmer), sizeof(m_last_kmer));
	m_flag = m_last_kmer >> (WLEN-2);
	m_last_kmer &= 0x3FFFFFFFFFFFFFFF; // The last two bits are flags

	m_kmers_count++;
}


void KmersSingleDataBaseSortedFile::load_kmers_upto_x(const uint64_t &threshold, std::vector<uint64_t> &kmers) {
	vector<uint64_t> temp_flags;
	load_kmers_upto_x(threshold, kmers, temp_flags);
}

void KmersSingleDataBaseSortedFile::load_kmers_upto_x(const uint64_t &threshold, std::vector<uint64_t> &kmers, std::vector<uint64_t> &flags) {
	if(!m_fin.is_open()) {
		throw std::logic_error("Trying to read k-mers from a non-open file.");
	}
	kmers.resize(0);
	flags.resize(0);
	while((m_last_kmer <= threshold) && (m_kmers_count<m_kmers_in_file)) {
		kmers.push_back(m_last_kmer);
		flags.push_back(m_flag);
		read_kmer(); // over-ride last kmer
	}
	if((m_last_kmer <= threshold) && (m_kmers_count==m_kmers_in_file)) {
		if(m_last_kmer != NULL_KEY) {
			flags.push_back(m_flag);
			kmers.push_back(m_last_kmer); // So we won't read the lasy k-mer twice
		}
		m_last_kmer = NULL_KEY;
		m_flag = NULL_KEY;
	}
}


