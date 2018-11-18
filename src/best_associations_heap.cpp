/**
 *       @file  best_associations_map.cpp
 *      @brief  
 *
 * Detailed description starts here.
 *
 *     @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
 *
 *   @internal
 *     Created  11/15/18
 *    Revision  $Id: doxygen.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
 *    Compiler  gcc/g++
 *     Company  Max Planck Institute for Developmental Biology Dep 6
 *   Copyright  Copyright (c) 2018, Yoav Voichek
 *
 * This source code is released for free distribution under the terms of the
 * GNU General Public License as published by the Free Software Foundation.
 * =====================================================================================
 */


#include "best_associations_heap.h"
using namespace std;
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
	kmer_score_priority_queue().swap(temp_queue); //make sure that temp is really emptied
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
	kmer_score_priority_queue().swap(temp_queue); //make sure that temp is really emptied
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
	kmer_score_priority_queue().swap(temp_queue); //make sure that temp is really emptied
	return res;
}

kmers_output_list kmer_heap::get_kmers_for_output(const size_t &kmer_len) const {
	kmers_output_list res;
	res.next_index = 0;

	kmer_score_priority_queue temp_queue(m_best_kmers);
	while(!temp_queue.empty()) {
		res.list.push_back(make_tuple(
					get<0>(temp_queue.top()),temp_queue.size(),
					get<2>(temp_queue.top())));
		temp_queue.pop();
	}
	kmer_score_priority_queue().swap(temp_queue); //make sure that temp is really emptied

	// sort list
	sort(begin(res.list), end(res.list), [](auto const &t1, auto const &t2) {
			return get<2>(t1) < get<2>(t2);});
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

