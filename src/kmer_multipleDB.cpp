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
#include <bitset>
#include <algorithm>
#include <numeric>
#include "../include/fisher-exact/kfunc.h"
using namespace std;

///
/// @brief  Ctor of kmer_multipleDB - 
/// @param  path_to_DBs:	a path in the filesystem to directories containing the files ok k-mers
//			db_names:		a list of names of all the DBs to use
//			sorted_kmer_fn:	filename inside each subdirectory containing the sorted k-mer list
/// @return 
///
kmer_multipleDB::kmer_multipleDB(
		const string& path_to_DBs, 
		const std::vector<std::string>& db_names, 
		const string& sorted_kmer_fn):
	m_DBs(),
	m_db_names(db_names),
	m_kmer_temp(), // just a temp vector
	m_accessions(db_names.size()), 
	m_kmers_pa(),
	m_hash_words((m_accessions+64-1)/64),
	m_verbose(true)
{
	// build all the kmer_DB objects and open the sorted k-mer file
	for(size_t i=0; i < m_db_names.size(); i++) {
		m_DBs.emplace_back(path_to_DBs +"/" + m_db_names[i], m_db_names[i]);
		m_DBs.back().open_sorted_kmer_file(sorted_kmer_fn);
	}
	if(m_verbose) 
		cerr << "We have " << m_accessions << " accessions and " << m_kmers_pa.size() << " hash-map\n";

	// Needs to define a "delete" key which is a non-possible input (our k-mer will be max 62 bits)
	m_kmers_pa.set_empty_key(0xFFFFFFFFFFFFFFFF);
}



///
/// @brief  filter from a kmer vector only the kmers part of a given set
/// @param  1. vector of uint64 representing k-mers and a set of kmers (hash set)
/// @return 
///
void filter_kmers_to_set(std::vector<uint64> &kmers, const kmer_set &set_kmers) {
	// need to implement....
	size_t move_to = 0;
	for(size_t i=1; i<kmers.size(); i++) {
		if(lookup_x(set_kmers,kmers[i])) {
			kmers[move_to] = kmers[i];
			move_to++;
		}
	}
	kmers.resize(move_to);
}


/**
 * @brief   kmer_multipleDB::load_kmers - load a subset of k-mers from sorted files
 * @param   Which iteration (iter) is it from all iterations (total_iter)
 *			as we go over sorted k-mer files, we iterate each time until a threshold
 *			Notice: the last 2-bits are 0 in all k-mers as k-mers are 31bp (62 bits)
 *			so the last threshold is 001111111...111 = 0x3FFFFF...FF
 * @return  
 */
void kmer_multipleDB::load_kmers(const uint64 &iter, const uint64 &total_iter, const kmer_set &set_kmers) {
	// as k-mers are 31 bp - the largest possible value is 0011111111...1111
	uint64 current_threshold = ((0x3FFFFFFFFFFFFFFF / total_iter)+1)*iter;
	m_kmers_pa.clear();
	my_multi_hash::iterator it_hash;
	for(size_t acc_i = 0; acc_i < m_accessions; ++acc_i) {
		size_t hashmap_i = acc_i / 64; // which word to use
		size_t bit_i = acc_i % 64; // which bit in word to use

		uint64 or_val = 1ull << bit_i; // create the word to use for modifying

		// Create the new word to put in new k-mers
		array<uint64, WORD64HASHT> new_bits;
		for(size_t i=0; i<WORD64HASHT; i++) new_bits[i] = 0ull;
		new_bits[hashmap_i] = or_val;

		// Reading new-kmers from file
		m_DBs[acc_i].read_sorted_kmers(m_kmer_temp, current_threshold); // m_kmer_temp is emptied in func'
		if(set_kmers.size() != 0) { 
			cerr << "using a kmer set! " << endl;
			cerr << "we have " << m_kmer_temp.size() << " kmers";
			filter_kmers_to_set(m_kmer_temp, set_kmers); // need to implement this
			cerr << " left with: " << m_kmer_temp.size() << endl;
		}

		// adding new k-mers info
		for(auto const& it: m_kmer_temp) {
			it_hash = m_kmers_pa.find(it);			
			if(it_hash == m_kmers_pa.end()) 
			{ m_kmers_pa.insert(my_multi_hash::value_type(it, new_bits)); 
			} else {it_hash->second[hashmap_i] |= or_val;}
		}
	}
}

uint64 reverseOne(uint64 x) {
	x = ((x & 0xFFFFFFFF00000000) >> 32) | ((x & 0x00000000FFFFFFFF) << 32);
	x = ((x & 0xFFFF0000FFFF0000) >> 16) | ((x & 0x0000FFFF0000FFFF) << 16);
	x = ((x & 0xFF00FF00FF00FF00) >> 8)  | ((x & 0x00FF00FF00FF00FF) << 8);
	x = ((x & 0xF0F0F0F0F0F0F0F0) >> 4)  | ((x & 0x0F0F0F0F0F0F0F0F) << 4);
	x = ((x & 0xCCCCCCCCCCCCCCCC) >> 2)  | ((x & 0x3333333333333333) << 2);
	x = ((x & 0xAAAAAAAAAAAAAAAA) >> 1)  | ((x & 0x5555555555555555) << 1);
	return x;
}

///
/// @brief  output all the kmers and the presence absence info found in class
//			output to stdout (can change to output to a given stream)
/// @param  
/// @return 
///
void kmer_multipleDB::output_kmers_textual() const { // output all k-mers found in the hash to stdout
	// last word length in bits
	size_t bits_last_word = m_accessions % 64;
	if(bits_last_word == 0) bits_last_word = 0;

	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) {		
		cout << bits2kmer31(it->first) << "\t"; 
		for(size_t i=0; i<(m_hash_words-1); i++) 
			cout << bitset<64>(reverseOne(it->second[i])).to_string();
		cout << bitset<64>(reverseOne(it->second[m_hash_words-1])).to_string() << endl;
		//			.substr(0,bits_last_word) << endl;
	}
}



///
/// @brief  output all the k-mers and presence/absence info in binary format
/// @param  filename to output the data to
/// @return 
///
void kmer_multipleDB::output_kmers_binary(const std::string &filename) const {
	ofstream of(filename, ios::binary);
	// First writing the file header
	of << (char)(0x31) << (char)(0xAD) << (char)(0x26) << (char)(0x00) ; // constant prefix
	of.write(reinterpret_cast<const char *>(&m_accessions), sizeof(m_accessions));			// size_t

	// writing the contents of the hash map
	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) {		
		of.write(reinterpret_cast<const char *>(&(it->first)), sizeof(it->first)); // write the k-mer	
		for(size_t i=0; i<m_hash_words; i++) 
			of.write(reinterpret_cast<const char *>(&(it->second[i])), sizeof(it->second[i])); 
	}
	of.close();
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
void kmer_multipleDB::output_plink_bed_file(const std::string &base_name) const  {
	ofstream f_bed(base_name + ".bed",ios::binary);
	ofstream f_bim(base_name + ".bim", ios::out);
	f_bed << (char)0x6C << (char)(0x1B) << (char)(0x01); // header of bed file

	unsigned char b;
	uint64 w;
	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) {		
		f_bim << "0\t" << bits2kmer31(it->first) << "\t0\t0\t0\t1\n"; 
		for(size_t i=0; i<m_hash_words; i++) {
			w = it->second[i];
			for(size_t bi=0; bi<16; bi++) { // every word is 64 accessions (16*4) every 4 accessions is a byte
				b = (w&1);
				w >>= 1;
				b ^= ((w&1)<<2);
				w >>= 1;
				b ^= ((w&1)<<4);
				w >>= 1;
				b ^= ((w&1)<<6);
				b |= (b<<1);
				w >>= 1;
				f_bed << b;
			}
		}
	}
	f_bed.close();
	f_bim.close();
}



///
/// @brief  go over all the kmers now present in the multipleDB, calculate association score and add this to
//			the given heap
/// @param  kmers_and_scores - kmers heap to save the best scores
//			scores - the measurments to associate the presence/absence with
//			names_scores - name of the DB the score is relevant to
/// @return 
///
// for continuous variable
void kmer_multipleDB::add_kmers_to_heap(kmer_heap &kmers_and_scores, const vector<double> &scores,
		const vector<string> &names_scores) const {
	vector<size_t> index_dbs = get_dbs_indices(names_scores); // finding the right indices of dbs
	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) 
		kmers_and_scores.add_kmer(it->first, calculate_kmer_score(it, scores, index_dbs));
}
// for categorical values (scores should have only 1 & 0!!!)
void kmer_multipleDB::add_kmers_to_heap(kmer_heap &kmers_and_scores, const vector<size_t> &scores,
		const vector<string> &names_scores) const {
	vector<size_t> index_dbs = get_dbs_indices(names_scores); // finding the right indices of dbs
	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) 
		kmers_and_scores.add_kmer(it->first, calculate_kmer_score(it, scores, index_dbs));
}


// return the indices of DB names inserted in the class DBs
// XX should we check if something is missing and raise exception?
vector<size_t> kmer_multipleDB::get_dbs_indices(const vector<string> &names) const {
	std::vector<size_t> indices ;
	for(auto n: names) 
		indices.push_back(std::find(m_db_names.begin(), m_db_names.end(), n) - m_db_names.begin());
	return indices;
}



///
/// @brief calculate the two sided t-test for association b/w phenotype (other score) to presence/absence
//			of a k-mer. 
/// @param  it - kmer and the presence/absence information
//			scores - vector of phenotype value to associate with
//			indices - the indices of the DBs the scores are for
//			min_in_group - give score only if at least this number is found in both cases
/// @return t-test statistics (or 0 if  min_in_group doesn't pass)
///
double kmer_multipleDB::calculate_kmer_score(
		const my_multi_hash::const_iterator& it, 
		const vector<double> &scores, 
		const vector<size_t> &indices,
		double min_in_group) const {
	double Ex0, Ex1, E2x0, E2x1, Sx0, Sx1, N0, N1, bit;
	Ex0 = Ex1 = Sx0 = Sx1 = N0 = N1 = E2x0 = E2x1 = 0;
	size_t i;
	for(size_t i_ind=0; i_ind<indices.size(); i_ind++) {
		i = indices[i_ind];
		size_t hashmap_i = i / 64;
		size_t bit_i = i % 64;
		bit = double((it->second[hashmap_i] >> bit_i)&1);
		Ex1 += bit*scores[i_ind];
		E2x1 += bit*scores[i_ind]*scores[i_ind];
		Ex0 += (1-bit)*scores[i_ind];
		E2x0 += (1-bit)*scores[i_ind]*scores[i_ind];
		N1 = N1+bit;
	}
	N0 = indices.size() - N1;
	if((min_in_group<=N0) && (min_in_group <= N1)) {
		Ex0 /= N0;
		E2x0 /= N0;
		Ex1 /= N1;
		E2x1 /= N1;

		Sx0 = sqrt(E2x0 - Ex0*Ex0);
		Sx1 = sqrt(E2x1 - Ex1*Ex1);	


		double v = N1 + N0 - 2; // degree of freedom
		double sp = sqrt(((N1-1) * Sx1 * Sx1 + (N0-1) * Sx0 * Sx0) / v); // Pooled variance
		return (Ex1 - Ex0) / (sp * sqrt(1.0 / N1 + 1.0 / N0)); // t-test statistics
	} else {return 0;}
}

///
/// @brief	calculate the fisher exat test for association b/w categorical phenotype (other score) to 
//			presence/absence of a k-mer. 
/// @param  it - kmer and the presence/absence information
//			scores - vector of phenotype value to associate with
//			indices - the indices of the DBs the scores are for
//			min_in_group - give score only if at least this number is found in both cases
/// @return ? statistics (or 0 if  min_in_group doesn't pass)
/// @note	fisher two-sided p-value this function can be more efficent if needed...
double kmer_multipleDB::calculate_kmer_score(
		const my_multi_hash::const_iterator& it,
		const vector<size_t> &scores, 
		const vector<size_t> &indices,
		double min_in_group) const {
	int n_all = indices.size();
	int n_presences = 0;
	int n_positive = accumulate(scores.begin(), scores.end(), 0); // ofcourse can be done outside...
	int n_presence_positive = 0;
	int bit;
	size_t i;
	for(size_t i_ind=0; i_ind<indices.size(); i_ind++) {
		i = indices[i_ind];
		size_t hashmap_i = i / 64;
		size_t bit_i = i % 64;
		bit = int((it->second[hashmap_i] >> bit_i)&1);
		n_presences += bit;
		n_presence_positive += bit*scores[i_ind];
	}
	if((n_presences <= min_in_group) || ((n_all-n_presences)<=min_in_group)) 
		return 0;
	int n11 = n_presence_positive;
	int n12 = n_positive - n_presence_positive;
	int n21 = n_presences - n_presence_positive;
	int n22 = n_all - (n11+n12+n21);
	double fisher_left_p, fisher_right_p, fisher_twosided_p;
	kt_fisher_exact(n11, n12, n21, n22,
			&fisher_left_p, &fisher_right_p, &fisher_twosided_p);
	//	cerr << "left greater - pval " << fisher_left_p << endl;
	//	cerr << "right greater - pval " << fisher_left_p << endl;
	//	cerr << "twosided - pval " << fisher_twosided_p << endl;
	return -log(fisher_twosided_p);
}


///
/// @brief  Ctor of kmer_heap initialized the priority queue
/// @param  what is the maximal number of results to save
/// @return 
///
kmer_heap::kmer_heap(size_t max_results):
	m_n_res(max_results),
	m_best_kmers(),
	cnt_kmers(0),
	cnt_pops(0) {
	}


///
/// @brief  add_kmer - add a new kmer to the heap
/// @param  
/// @return 
///
void kmer_heap::add_kmer(const uint64 &k, const double &score) {
	const static cmp_second compare_func;
	cnt_kmers++;
	if(m_best_kmers.size() < m_n_res) {
		m_best_kmers.push(kmer_score(k, score));
	} else {
		kmer_score new_res(k, score);
		if(compare_func(new_res, m_best_kmers.top())) {
			cnt_pops++;
			m_best_kmers.pop();
			m_best_kmers.push(new_res);
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
		of.write(reinterpret_cast<const char *>(&temp_queue.top().first), sizeof(uint64));
		temp_queue.pop();
	}
	of.close();
}


///
/// @brief  output the k-mers with the scores to a file
/// @param  
/// @return 
///
void kmer_heap::output_to_file_with_scores(const std::string &filename) const {
	kmer_score_priority_queue temp_queue(m_best_kmers);
	ofstream of(filename, ios::binary);
	while(!temp_queue.empty()) {
		of.write(reinterpret_cast<const char *>(&temp_queue.top().first), sizeof(uint64));
		of.write(reinterpret_cast<const char *>(&temp_queue.top().second), sizeof(double));
		temp_queue.pop();
	}
	of.close();
}
