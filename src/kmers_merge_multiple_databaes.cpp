/**
 *       @file  kmers_merge_multiple_databaes.cpp
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

#include "kmers_merge_multiple_databaes.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <algorithm>
#include <bitset>
#include <numeric>

using namespace std;


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

