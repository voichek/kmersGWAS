///
///      @file  kmer_multipleDB.cpp
///     @brief  
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  07/19/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmer_multipleDB.h"

#include <math.h>
#include "nmmintrin.h" //_mm_popcnt_u64 
#include <algorithm>
#include <bitset>
#include <numeric>
using namespace std;

///
/// @brief  Ctor of kmer_multipleDB - 
/// @param 	DB_paths:	paths in the filesystem to directories containing the files of k-mers
//		db_names:	a list of names of all the DBs to use (same order as DB_paths)
//		sorted_kmer_fn:	filename inside each subdirectory containing the sorted k-mer list
/// @return 
///
kmer_multipleDB::kmer_multipleDB(
		const string &merge_db_file,
		const vector<string> &db_names,
		const vector<string> &db_to_use, 	
		const uint32 &kmer_len):
	m_db_names_db_file(db_names), // table in file
	m_db_names_table(db_to_use), // table in memory
	m_accessions_db_file(db_names.size()), 
	m_accessions(db_to_use.size()),
	m_kmer_table_file(merge_db_file, ios::binary | ios::ate), // ios::ate - put pointer in the end of file
	m_hash_words_db_file((m_accessions_db_file+WLEN-1)/WLEN),
	m_hash_words((m_accessions+WLEN-1)/WLEN),
	m_kmers(),
	m_kmers_table(),
	m_kmer_len(kmer_len),
	m_left_in_file(),
	m_kmer_number(),
	m_kmer_loaded(0),
	m_row_offset(0),
	m_map_word_index(),
	m_map_bit_index(),
	m_verbose(true)
{
	if(m_kmer_table_file.is_open()) {
		m_left_in_file = m_kmer_table_file.tellg(); // get file size
		if(m_left_in_file <= (4 + 8 + 4)) 
			throw std::logic_error("Kmer table size is too small");
		m_kmer_table_file.seekg (0, ios::beg); //go back to begining

		uint32 prefix, file_kmer_len;
		uint64_t file_accession_number;
		m_kmer_table_file.read(reinterpret_cast<char *>(&prefix), sizeof(prefix));
		m_kmer_table_file.read(reinterpret_cast<char *>(&file_accession_number), sizeof(file_accession_number));
		m_kmer_table_file.read(reinterpret_cast<char *>(&file_kmer_len), sizeof(file_kmer_len));
		m_left_in_file -= (sizeof(prefix) + sizeof(file_accession_number) + sizeof(file_kmer_len));

		if(prefix!=0xDDCCBBAA)
			throw std::logic_error("Incorrect prefix");
		if(file_accession_number != m_accessions_db_file)
			throw std::logic_error("Number of accession in file not as defined in class");
		if(file_kmer_len != m_kmer_len)
			throw std::logic_error("Kmer length not as defined in class");

		// From the size of the file we can calculate the number of k_mers
		size_t size_per_kmer = sizeof(uint64_t) * (1 + m_hash_words_db_file);
		if((m_left_in_file % size_per_kmer) != 0)
			throw std::logic_error("size of file not valid");
		m_kmer_number = m_left_in_file / size_per_kmer;
		cerr << "We have " << m_kmer_number << endl;
		create_map_from_all_DBs();
	} else {
		throw std::logic_error("Couldn't open kmer table file: " + merge_db_file);
	}
}



/**
 * @brief   Load k-mers from table to memory
 * The function will also squeeze to keep only k-mers which are part of the current DB set		
 * @param   size of memory batch and a set of k-mers to intersect with
 * @return return false if table file is allready finished 
 */
bool kmer_multipleDB::load_kmers(const uint64_t &batch_size, const kmer_set &set_kmers_to_use) {
	m_row_offset = m_kmer_loaded;
	clear();
	if(m_left_in_file == 0) {
		return false;
		m_kmer_table_file.close();
	} else { // file not empty
		vector<uint64_t> reading_buffer(m_hash_words_db_file+1);
		for(size_t i=0; (i<batch_size) && (m_left_in_file>0); i++) {
			// Reading from file
			//		for(size_t j=0;j<reading_buffer.size();j++)
			//			m_kmer_table_file.read(reinterpret_cast<char *>(&reading_buffer[j]), 
			//					sizeof(reading_buffer[j]));
			m_kmer_table_file.read(reinterpret_cast<char *>(reading_buffer.data()), 
					sizeof(uint64_t)*reading_buffer.size());
			//update counters
			m_left_in_file -= (reading_buffer.size() * sizeof(uint64_t));
			m_kmer_loaded++;
			if((set_kmers_to_use.size() == 0) || (lookup_x(set_kmers_to_use,reading_buffer[0]))) {
				m_kmers.push_back(reading_buffer[0]); // save k-mers

				// add info to table in memory (= squeeze)
				size_t table_offset = m_kmers_table.size();
				m_kmers_table.resize(table_offset+m_hash_words, 0);
				for(size_t col_index=0; col_index<m_accessions; col_index++) {
					uint64_t new_bit = (reading_buffer[m_map_word_index[col_index]+1] >> m_map_bit_index[col_index]) & 1ull;
					uint64_t hashmap_i = (col_index>>6);
					uint64_t  bit_i = col_index&63;
					m_kmers_table[table_offset + hashmap_i] |= (new_bit << bit_i);
				}
			}	
		}
	}
	return true;
}


///
/// @brief  output all the kmers and the presence absence info found in class
//			output to stdout (can change to output to a given stream)
/// @param  
/// @return 
///
void kmer_multipleDB::output_kmers_textual() const { // output all k-mers found in the hash to stdout
	// last word length in bits
	for(size_t kmer_i=0; kmer_i<m_kmers.size(); kmer_i++) {
		cout << bits2kmer31(m_kmers[kmer_i], m_kmer_len) << "\t";
		size_t container_i = kmer_i*m_hash_words;
		for(size_t i=container_i; i<(container_i+m_hash_words-1); i++) 
			cout << bitset<WLEN>(reverseOne(m_kmers_table[i])).to_string();
		cout << bitset<WLEN>(reverseOne(m_kmers_table[container_i+m_hash_words-1])).to_string() << endl;
	}
}

///
/// @brief  output all k-mers and presence/absence information in a format readable by Plink
//			that is "bed" format. We will also create a meta file of the k-mers in the bed file.
//
//			Format taken from the plink website filetypes section:
//			https://www.cog-genomics.org/plink/1.9/formats#bed
//			"The first three bytes should be 0x6c, 0x1b, and 0x01 in that order. (There are old versions of 
//			the .bed format which start with a different "magic number"; PLINK 1.9 recognizes them, but will 
//			convert sample-major files to the current variant-major format on sight. See the bottom of the 
//			original .bed definition page for details; that page also contains a more verbose version of the 
//			discussion below.)
//
//			The rest of the file is a sequence of V blocks of N/4 (rounded up) bytes each, where V is the 
//			number of variants and N is the number of samples. The first block corresponds to the first marker 
//			in the .bim file, etc.
//
//			The low-order two bits of a block's first byte store the first sample's genotype code. ("First 
//			sample" here means the first sample listed in the accompanying .fam file.) The next two bits 
//			store the second sample's genotype code, and so on for the 3rd and 4th samples. The second byte 
//			stores genotype codes for the 5th-8th samples, the third byte stores codes for the 9th-12th, etc.
//
//			The two-bit genotype codes have the following meanings:
//
//			00    Homozygous for first allele in .bim file
//			01    Missing genotype
//			10    Heterozygous
//			11    Homozygous for second allele in .bim file
//
/// @param  base path to write the .bin & .bed file
/// @return 
/// @notes	I might need to add the logic of taking out duplicates patterns here...
void kmer_multipleDB::output_plink_bed_file(const string &base_name) const  {
	bedbim_handle f_handle(base_name);
	output_plink_bed_file(f_handle);
	f_handle.close();
};

void kmer_multipleDB::output_plink_bed_file(bedbim_handle &f, const kmer_set &set_kmers) const  {
	for(size_t kmer_i=0; kmer_i<m_kmers.size(); kmer_i++) {
		if((set_kmers.size() == 0) || (lookup_x(set_kmers,m_kmers[kmer_i]))) { // check k-mer in set (or empty set)
			write_PA(bits2kmer31(m_kmers[kmer_i], m_kmer_len), kmer_i, f); 
		}
	}
}

void kmer_multipleDB::write_PA(const string &name, const size_t &kmer_i, bedbim_handle &f) const {
	f.f_bim << "0\t" << name << "\t0\t0\t0\t1\n"; 
	size_t container_i = kmer_i*m_hash_words;
	size_t acc_index = 0;
	for(size_t i=container_i; i<(container_i+m_hash_words); i++) {
		uint64_t w = m_kmers_table[i];
		for(size_t bi=0; (bi<16) && (acc_index < m_accessions); bi++) { // every word is 64 accessions (16*4) every 4 accessions is a byte
			unsigned char b = (w&1);
			w >>= 1;
			b ^= ((w&1)<<2);
			w >>= 1;
			b ^= ((w&1)<<4);
			w >>= 1;
			b ^= ((w&1)<<6);
			b |= (b<<1);
			w >>= 1;
			f.f_bed << b;
			acc_index += 4;
		}
	}

}

size_t kmer_multipleDB::output_plink_bed_file(bedbim_handle &f, const vector<kmer_output> &kmer_list, size_t index) const {
	for(size_t kmer_i=0; kmer_i<m_kmers.size(); kmer_i++) {
		if((index<kmer_list.size()) && 
				(get<1>(kmer_list[index])==(kmer_i+m_row_offset))) {
			write_PA(get<0>(kmer_list[index]), kmer_i, f);
			index++;	
		}
	}
	return index;
}

/////
///// @brief  go over all the kmers now present in the multipleDB, calculate association score and add this to
////			the given heap
///// @param  kmers_and_scores - kmers heap to save the best scores
////			scores - the measurments to associate the presence/absence with
////			names_scores - name of the DB the score is relevant to
///// @return 
/////
void kmer_multipleDB::add_kmers_to_heap(kmer_heap &kmers_and_scores, vector<float> scores, 
		const size_t &min_cnt) const {
	// Due to efficency consideration we pass scores not by reference and change it to be a multiplication of
	// word size (saving many not neccessery "if" in the scoring procedure)
	scores.resize(m_hash_words*sizeof(uint64_t)*8, 0); 
	float sum_scores(0);
	for(size_t i=0; i<scores.size(); i++) 
		sum_scores += scores[i];
	for(size_t kmer_index=0; kmer_index<m_kmers.size(); kmer_index++)
		kmers_and_scores.add_kmer(m_kmers[kmer_index], 
				calculate_kmer_score(kmer_index, scores, sum_scores,  min_cnt),
				m_row_offset + kmer_index);
}


// return the indices of DB names inserted in the class DBs
void kmer_multipleDB::create_map_from_all_DBs() {
	m_map_word_index.resize(0);
	m_map_bit_index.resize(0);
	size_t index_file_table;
	for(size_t index_mem_table=0; index_mem_table < m_accessions; index_mem_table++) {
		index_file_table = find(m_db_names_db_file.begin(), m_db_names_db_file.end(), m_db_names_table[index_mem_table]) 
			- m_db_names_db_file.begin();
		if(index_file_table == m_accessions_db_file) // Couldn't find it
			throw std::logic_error("All accessions suppose to be in DB file: " + m_db_names_table[index_mem_table]);
		m_map_word_index.push_back(index_file_table / WLEN);
		m_map_bit_index.push_back(index_file_table % WLEN);
	}
}

#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
double kmer_multipleDB::calculate_kmer_score(
		const size_t kmer_index, 
		const vector<float> &scores, // Scores has to be a multiplication of wordsize (64) 
		const float score_sum,
		const uint64_t min_in_group
		) const {
	uint64_t bit, N1(0), N0;
	float Ex1(0);
	double N=(double)m_accessions;
	size_t container_i = kmer_index*m_hash_words;
	/* The following section is where my program is most of the time */
	size_t i=0;
	for(uint64_t hashmap_i=0; hashmap_i<m_hash_words; hashmap_i++) {
		uint64_t word = m_kmers_table[container_i+hashmap_i];
		N1 += _mm_popcnt_u64(word);
		for(uint64_t bit_i=0; bit_i<(8*sizeof(uint64_t)); bit_i++) {
			bit = word & 1;
			Ex1 += (scores[i]*bit);
			word = (word>>1);		
			i++;
		}
	}
	/* End of critical section */
	N0 = m_accessions - N1;
	if((min_in_group<=N0) && (min_in_group <= N1)) {
		double sum_gi = (double)N1;
		double yigi = (double)Ex1;
		double r =  N*yigi- sum_gi*score_sum;
		r = r * r;
		return r / (N*sum_gi - sum_gi*sum_gi);
	} else {return 0;}
}
#pragma GCC pop_options

///
/// @brief  Ctor of kmer_heap initialized the priority queue
/// @param  what is the maximal number of results to save
/// @return 
///
kmer_heap::kmer_heap(size_t max_results):
	m_n_res(max_results),
	m_best_kmers(),
	cnt_kmers(0),
	cnt_pops(0),
	cnt_push(0),
	lowest_score(0){
	}


///
/// @brief  add_kmer - add a new kmer to the heap
/// @param  k - kmer
//			score - kmers score
//			kmer_row - number of the row in kmers table (will be useful for retrieving the kmer)
/// @return 
///
void kmer_heap::add_kmer(const uint64_t &k, const double &score, const uint64_t &kmer_row) {
	cnt_kmers++;
	if(m_best_kmers.size() < m_n_res) { // if heap is not full
		m_best_kmers.push(kmer_score(k, score, kmer_row));
		cnt_push++;
		lowest_score = get<1>(m_best_kmers.top());
	} else {
		if(score > lowest_score) {
			kmer_score new_res(k, score, kmer_row);
			cnt_pops++;
			cnt_push++;
			m_best_kmers.pop();
			m_best_kmers.push(new_res);
			lowest_score = get<1>(m_best_kmers.top());
		}
	}
}


///
/// @brief  output only the k-mers to a file (for future use)
/// @param  
/// @return 
///
void kmer_heap::output_to_file(const string &filename) const {
	kmer_score_priority_queue temp_queue(m_best_kmers);
	ofstream of(filename, ios::binary);
	while(!temp_queue.empty()) {
		of.write(reinterpret_cast<const char *>(&get<0>(temp_queue.top())), sizeof(uint64_t));
		temp_queue.pop();
	}
	of.close();
}


///
/// @brief  output the k-mers with the scores to a file
///
void kmer_heap::output_to_file_with_scores(const std::string &filename) const {
	kmer_score_priority_queue temp_queue(m_best_kmers);
	ofstream of(filename, ios::binary);
	while(!temp_queue.empty()) {
		of.write(reinterpret_cast<const char *>(&get<0>(temp_queue.top())), sizeof(uint64_t));
		of.write(reinterpret_cast<const char *>(&get<1>(temp_queue.top())), sizeof(double));
		temp_queue.pop();
	}
	of.close();
}


/// @brief  output all the k-mers in the heap
/// @return a kmer_set with all the k-mers from the heap
kmer_set kmer_heap::get_kmer_set() const {
	kmer_set res;
	res.set_empty_key(NULL_KEY); // need to define empty value for google dense hash table

	kmer_score_priority_queue temp_queue(m_best_kmers);
	while(!temp_queue.empty()) {
		res.insert(get<0>(temp_queue.top()));
		temp_queue.pop();
	}
	return res;
}
kmers_output_list kmer_heap::get_kmers_for_output(const size_t &kmer_len) const {
	kmers_output_list res;
	res.next_index = 0;
	
	kmer_score_priority_queue temp_queue(m_best_kmers);
	while(!temp_queue.empty()) {
		res.list.push_back(make_tuple(
					bits2kmer31(get<0>(temp_queue.top()),kmer_len) + "_" + to_string(temp_queue.size()),
					get<2>(temp_queue.top())));
		temp_queue.pop();
	}
	
	// sort list
	sort(begin(res.list), end(res.list), [](auto const &t1, auto const &t2) {
			        return get<1>(t1) < get<1>(t2);});
	return res;
}	

/// @bried	output the status (pop/ push/ size) of the heap to stderr
void kmer_heap::plot_stat() const { 
	std::cerr << "[heap-stat] max\t" << m_n_res <<
		"\tsize\t" << m_best_kmers.size() <<
		"\tkmers\t" << cnt_kmers <<
		"\tpops\t" << cnt_pops <<
		"\tpush\t" << cnt_push << endl; 
}

// Close handles of bed & bim files
void bedbim_handle::close() {
	if(f_bed.is_open())
		f_bed.close();
	if(f_bim.is_open())
		f_bim.close();
}


kmer_multipleDB_merger::kmer_multipleDB_merger(const vector<string> &DB_paths,
		const vector<string> &db_names, 
		const string &sorted_kmer_fn,
		const uint32 &kmer_len):
	m_DBs(),
	m_db_names(db_names),
	m_kmer_temp(), // just a temp vector
	m_accessions(db_names.size()), 
	m_hash_words((m_accessions+WLEN-1)/WLEN),
	kmers_to_index(),
	container(),
	m_kmer_len(kmer_len)
{
	// build all the kmer_DB objects and open the sorted k-mer file
	for(size_t i=0; i < m_db_names.size(); i++) {
		cerr << "Open DB: " << m_db_names[i] << endl;
		m_DBs.emplace_back(DB_paths[i] , m_db_names[i], m_kmer_len);
		m_DBs.back().open_sorted_kmer_file(sorted_kmer_fn);
	}

	// Needs to define a "delete" key which is a non-possible input (our k-mer will be max 62 bits)
	kmers_to_index.set_empty_key(NULL_KEY);
}

void kmer_multipleDB_merger::output_table_header(ofstream& T) const {
	T << (char)(0xAA) << (char)(0xBB) << (char)(0xCC) << (char)(0xDD) ; // constant prefix
	T.write(reinterpret_cast<const char *>(&m_accessions), sizeof(m_accessions));			// size_t
	T.write(reinterpret_cast<const char *>(&m_kmer_len),   sizeof(m_kmer_len));
}

void kmer_multipleDB_merger::output_to_table(ofstream& T) const {
	// writing the contents of the hash map
	for(my_hash::const_iterator it=kmers_to_index.begin(); it != kmers_to_index.end(); ++it) {		
		T.write(reinterpret_cast<const char *>(&(it->first)), sizeof(it->first)); // write the k-mer	
		for(size_t i=it->second; i<(it->second+m_hash_words); i++) 
			T.write(reinterpret_cast<const char *>(&(container[i])), sizeof(container[i])); 
	}
}

void kmer_multipleDB_merger::clear_content() {
	kmers_to_index.clear();
	container.resize(0);
} // clear container 


/**
 * @brief   kmer_multipleDB::load_kmers - load a subset of k-mers from sorted files
 * @param   Which iteration (iter) is it from all iterations (total_iter)
 *			as we go over sorted k-mer files, we iterate each time until a threshold
 *			Notice: the last 2-bits are 0 in all k-mers as k-mers are 31bp (62 bits)
 *			so the last threshold is 001111111...111 = 0x3FFFFF...FF
 * @return  
 */
void kmer_multipleDB_merger::load_kmers(const uint64_t &iter, const uint64_t &total_iter, const kmer_set &set_kmers) {
	// as k-mers are 31 bp - the largest possible value is 0011111111...1111
	uint64_t current_threshold = ((1ull << (m_kmer_len*2ull))-1ull);
	current_threshold = ((current_threshold / total_iter)+1)*iter;
	cerr << iter << " / " << total_iter << "\t:\t" << bitset<WLEN>(current_threshold) << endl;

	clear_content();
	my_hash::iterator it_hash;

	uint64_t index; 
	for(size_t acc_i = 0; acc_i < m_accessions; ++acc_i) {
		cerr << total_iter << "\t" << iter << "\t" << acc_i << endl; 
		size_t hashmap_i = acc_i / WLEN; // which word to use
		size_t bit_i = acc_i % WLEN; // which bit in word to use

		uint64_t or_val = 1ull << bit_i; // create the word to use for modifying

		// Reading new-kmers from file
		m_DBs[acc_i].read_sorted_kmers(m_kmer_temp, current_threshold); // m_kmer_temp is emptied in func'
		if(set_kmers.size() != 0)  
			filter_kmers_to_set(m_kmer_temp, set_kmers);

		// adding new k-mers info
		for(auto const& it: m_kmer_temp) {
			it_hash = kmers_to_index.find(it); // find if k_mer is allready in the hash map			
			if(it_hash == kmers_to_index.end()) { //if not
				kmers_to_index.insert(my_hash::value_type(it, container.size()));
				index = container.size() + hashmap_i;
				for(size_t i=0; i<m_hash_words; i++)
					container.push_back(0);
			} else {
				index = (it_hash->second)+hashmap_i;
			}
			container[index] |= or_val; //update container
		}
	}
}

