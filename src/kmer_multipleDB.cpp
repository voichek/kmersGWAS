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
/// @param 	DB_paths:	paths in the filesystem to directories containing the files of k-mers
//		db_names:	a list of names of all the DBs to use (same order as DB_paths)
//		sorted_kmer_fn:	filename inside each subdirectory containing the sorted k-mer list
/// @return 
///
kmer_multipleDB::kmer_multipleDB(
		const vector<string> &DB_paths, 
		const vector<string> &db_names, 
		const string &sorted_kmer_fn):
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
		cerr << m_db_names[i] << endl;
		m_DBs.emplace_back(DB_paths[i] , m_db_names[i]);
		m_DBs.back().open_sorted_kmer_file(sorted_kmer_fn);
	}

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
	cerr << iter << " / " << total_iter << "\t:\t" << bitset<64>(current_threshold) << endl;
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
			filter_kmers_to_set(m_kmer_temp, set_kmers); // need to implement this
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
void kmer_multipleDB::output_plink_bed_file(const string &base_name) const  {
	bedbim_handle f_handle(base_name);
	output_plink_bed_file(f_handle);
	f_handle.close();
};

void kmer_multipleDB::output_plink_bed_file(bedbim_handle &f, const kmer_set &set_kmers) const  {
	unsigned char b;
	uint64 w;
	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) {
		if((set_kmers.size() == 0) || (lookup_x(set_kmers,it->first))) { // check k-mer in set (or empty set)
			f.f_bim << "0\t" << bits2kmer31(it->first) << "\t0\t0\t0\t1\n"; 
			size_t acc_index = 0;
			for(size_t i=0; i<m_hash_words; i++) {
				w = it->second[i];
				for(size_t bi=0; (bi<16) && (acc_index < m_accessions); bi++) { // every word is 64 accessions (16*4) every 4 accessions is a byte
					b = (w&1);
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
	}
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

// efficent version;
void kmer_multipleDB::add_kmers_to_heap(kmer_heap &kmers_and_scores, const vector<double> &scores) const {	
	double sum_scores, sum_scores2;
	sum_scores = sum_scores2 = 0;
	vector<double> scores2(scores);
	for(size_t i=0; i<scores.size(); i++) {
		scores2[i] = scores2[i]*scores2[i];
		sum_scores += scores[i];
		sum_scores2 += scores2[i];
	}
	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) 
		kmers_and_scores.add_kmer(it->first, calculate_kmer_score(it, scores, scores2, sum_scores, sum_scores2));
}
//// for categorical values (scores should have only 1 & 0!!!)
//void kmer_multipleDB::add_kmers_to_heap(kmer_heap &kmers_and_scores, const vector<uint64> &scores, 
//		const size_t &min_cnt) const {
//	cerr << "[XXX] " << min_cnt << endl;
//	uint64 sum_scores, sum_scores2;
//	sum_scores = sum_scores2 = 0;
//	vector<uint64> scores2(scores);
//	for(size_t i=0; i<scores.size(); i++) {
//		scores2[i] = scores2[i]*scores2[i];
//		sum_scores += scores[i];
//		sum_scores2 += scores2[i];
//	}
//	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) 
//		kmers_and_scores.add_kmer(it->first, calculate_kmer_score(it, scores, scores2, sum_scores, sum_scores2, min_cnt));
//}

// 
void kmer_multipleDB::add_kmers_to_heap(kmer_heap &kmers_and_scores, const vector<uint64> &scores, 
		const size_t &min_cnt) const {
	uint64 sum_scores(0);
	for(size_t i=0; i<scores.size(); i++) 
		sum_scores += scores[i];
	double d_sum_scores = (double)sum_scores;
	for(my_multi_hash::const_iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) 
		kmers_and_scores.add_kmer(it->first, calculate_kmer_score(it, scores, d_sum_scores,  min_cnt));
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
		const double min_in_group) const {
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

// Precaculated scores^2 and sum of scores - imporve in efficency
double kmer_multipleDB::calculate_kmer_score(
		const my_multi_hash::const_iterator& it, 
		const vector<double> &scores, 
		const vector<double> &scores2,
		const double score_sum,
		const double score2_sum,
		const double min_in_group
		) const {
	double Ex0, Ex1, E2x0, E2x1, Sx0, Sx1, N0, N1, bit;
	Ex0 = Ex1 = Sx0 = Sx1 = N0 = N1 = E2x0 = E2x1 = 0;
	for(size_t i=0; i<scores.size(); i++) {
		size_t hashmap_i = i / 64;
		size_t bit_i = i % 64;
		bit = double((it->second[hashmap_i] >> bit_i)&1);
		Ex1 += bit*scores[i];
		E2x1 += bit*scores2[i];
		N1 = N1+bit;
	}
	Ex0 = score_sum - Ex1;
	E2x0 = score2_sum-E2x1;
	N0 = scores.size() - N1;
	if((min_in_group<=N0) && (min_in_group <= N1)) {
		Ex0 /= N0;
		E2x0 /= N0;
		Ex1 /= N1;
		E2x1 /= N1;
		double d = (Ex1-Ex0);
		d = d*d;

		double S2x0 = E2x0 - Ex0*Ex0;
		double S2x1 = E2x1 - Ex1*Ex1;	

		double sp2 = ((N1-1) * S2x1  + (N0-1) * S2x0) ; // Pooled variance
		return d / (sp2 * (1.0 / N1 + 1.0 / N0)); // t-test statistics

	} else {return 0;}
}

// 2nd efficency improvments - calculations in uints (+- are more efficient than doubles, but */ is less
// in total it does improve efficency
double kmer_multipleDB::calculate_kmer_score(
		const my_multi_hash::const_iterator& it, 
		const vector<uint64> &scores, 
		const vector<uint64> &scores2,
		const uint64 score_sum,
		const uint64 score2_sum,
		const uint64 min_in_group
		) const {
	uint64 Ex0, Ex1, E2x0, E2x1, N0, N1, bit;
	Ex0 = Ex1 =  N0 = N1 = E2x0 = E2x1 = 0;
	for(uint64 i=0; i<scores.size(); i++) {
		uint64 hashmap_i = (i>>6);
		uint64  bit_i = i&63;

		bit = (it->second[hashmap_i] >> bit_i)&1;
		N1 = N1+bit;
		bit = -bit; // so bit will be a mask

		Ex1 += scores[i]&bit;
		E2x1 += scores2[i]&bit;
	}
	Ex0 =  score_sum - Ex1;
	E2x0 = score2_sum - E2x1;
	N0 = scores.size() - N1;
	if((min_in_group<=N0) && (min_in_group <= N1)) {
		double dN0 = (double)N0;
		double dN1 = (double)N1;

		double dEx0  = (double)Ex0  / dN0;
		double dE2x0 = (double)E2x0 / dN0;
		double dEx1  = (double)Ex1  / dN1;
		double dE2x1 = (double)E2x1 / dN1;
		double d = (dEx1-dEx0);
		d = d*d;
		double S2x0 = dE2x0 - dEx0*dEx0;
		double S2x1 = dE2x1 - dEx1*dEx1;	

		double sp2 = ((dN1-1) * S2x1  + (dN0-1) * S2x0) ; // Pooled variance
		return d / (sp2 * (1.0 / dN1 + 1.0 / dN0)); // t-test statistics
	
	} else {return 0;}
}

double kmer_multipleDB::calculate_kmer_score(
		const my_multi_hash::const_iterator& it, 
		const vector<uint64> &scores, 
		const double score_sum,
		const uint64 min_in_group
		) const {
	uint64 bit, N1(0), Ex1(0), N0;
	double N=(double)scores.size();
	for(uint64 i=0; i<scores.size(); i++) {
		uint64 hashmap_i = (i>>6);
		uint64  bit_i = i&63;

		bit = (it->second[hashmap_i] >> bit_i)&1;
		N1 = N1+bit;
		bit = -bit; // so bit will be a mask

		Ex1 += scores[i]&bit;
	}
	N0 = scores.size() - N1;
	if((min_in_group<=N0) && (min_in_group <= N1)) {
		double sum_gi = (double)N1;
		double yigi = (double)Ex1;
		double r =  N*yigi- sum_gi*score_sum;
		r = r * r;
		return r / (N*sum_gi - sum_gi*sum_gi);
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
	cnt_pops(0),
	cnt_push(0) {
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
		cnt_push++;
	} else {
		kmer_score new_res(k, score);
		if(compare_func(new_res, m_best_kmers.top())) {
			cnt_pops++;
			cnt_push++;
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


/// @brief  output all the k-mers in the heap
/// @return a kmer_set with all the k-mers from the heap
kmer_set kmer_heap::get_kmer_set() const {
	kmer_set res;
	res.set_empty_key(NULL_KEY); // need to define empty value for google dense hash table

	kmer_score_priority_queue temp_queue(m_best_kmers);
	while(!temp_queue.empty()) {
		res.insert(temp_queue.top().first);
		temp_queue.pop();
	}
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
