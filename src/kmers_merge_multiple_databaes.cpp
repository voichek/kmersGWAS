///
///      @file  kmers_merge_multiple_databaes.cpp
///     @brief  Implementation of MultipleKmersDataBasesMerger 
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  11/15/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
///This source code is released for free distribution under the terms of the
///GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmers_merge_multiple_databaes.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <algorithm>
#include <bitset>
#include <numeric>

using namespace std;

MultipleKmersDataBasesMerger::MultipleKmersDataBasesMerger(
		const vector<string> &sorted_kmers_filenames,
		const vector<string> &accessions_names, 
		const string &all_possible_kmers_filename,
		const uint32 &kmer_len):
	m_sorted_kmers_list(),
	m_accessions_name(accessions_names),
	m_kmer_temp(), // just a temp vector
	m_accessions(accessions_names.size()), 
	m_hash_words((m_accessions+WLEN-1)/WLEN),
	m_kmers_to_index(),
	m_container(),
	m_container_kmers(),
	m_kmer_len(kmer_len),
	m_possible_kmers(all_possible_kmers_filename)
{
	// build all the KmersSingleDataBase objects and open the sorted k-mer file
	for(size_t i=0; i < m_accessions_name.size(); i++)
		m_sorted_kmers_list.emplace_back(sorted_kmers_filenames[i]);

	// Needs to define a "delete" key which is a non-possible input (our k-mer will be max 62 bits)
	m_kmers_to_index.set_empty_key(NULL_KEY);
}

void MultipleKmersDataBasesMerger::output_table_header(ofstream& T) const {
	T << (char)(0xAA) << (char)(0xBB) << (char)(0xCC) << (char)(0xDD) ; // constant prefix
	T.write(reinterpret_cast<const char *>(&m_accessions), sizeof(m_accessions));			// size_t
	T.write(reinterpret_cast<const char *>(&m_kmer_len),   sizeof(m_kmer_len));
}

void MultipleKmersDataBasesMerger::output_to_table(ofstream& T) const {
	// writing the contents of the hash map
	size_t index = 0;
	for(size_t kmer_index=0; kmer_index<m_container_kmers.size(); kmer_index++) {
		T.write(reinterpret_cast<const char *>(&(m_container_kmers[kmer_index])), 
				sizeof(m_container_kmers[kmer_index])); // write the k-mer
		for(size_t i=0; i<m_hash_words; i++) {
			T.write(reinterpret_cast<const char *>(&(m_container[index])), sizeof(m_container[index]));
			index++;
		}
	}
	cerr << "Wrote: kmers=" << m_container_kmers.size() << "\tpa words=" << index << "\tcontainer size=" << m_container.size()
		<< "\thash-map size=" << m_kmers_to_index.size() << endl;
}



/**
 * @brief   load a subset of k-mers from sorted files and code presence/absence in table
 * @param   Which iteration (iter) is it from all iterations (total_iter)
 *	    as we go over sorted k-mer files, we iterate each time until a threshold
 *	    Note 1: the last 2-bits are 0 in all k-mers as k-mers are 31bp (62 bits)
 *	    so the last threshold is 001111111...111 = 0x3FFFFF...FF
 *	    Note 2: We use hash-table to find relevant row in the table, as the two
 *	    lists are sorted this is not the most efficent solution.
 *	    (Might change in the future)
 * @return  
 */
void MultipleKmersDataBasesMerger::load_kmers(const uint64_t &iter, const uint64_t &total_iter) { 
	// clear relevant containers
	m_container.resize(0); 
	m_kmers_to_index.clear();

	// define threshold to read kmers
	uint64_t current_threshold = kmers_step_to_threshold(iter, total_iter, m_kmer_len);
	cerr << iter << " / " << total_iter << "\t:\t" << bitset<WLEN>(current_threshold) << endl;
	
	// First reads all the kmers from the general list
	m_possible_kmers.load_kmers_upto_x(current_threshold, m_container_kmers);
	m_container.resize(m_hash_words * m_container_kmers.size(), 0);
	
	// Intialize dictionary from kmers to index of container
	for(size_t kmer_index=0; kmer_index<m_container_kmers.size(); kmer_index++) 
		m_kmers_to_index.insert(KmerUint64Hash::value_type(m_container_kmers[kmer_index], kmer_index*m_hash_words));

	KmerUint64Hash::iterator it_hash;
	for(size_t acc_i = 0; acc_i < m_accessions; ++acc_i) { // Go over all accessions
		size_t hashmap_i = acc_i / WLEN; // which word to use
		size_t bit_i = acc_i % WLEN; 	 // which bit in word to use
		uint64_t or_val = 1ull << bit_i; // create the word to use for modifying

		// Reading kmers from file
		m_sorted_kmers_list[acc_i].load_kmers_upto_x(current_threshold, m_kmer_temp);

		// adding new k-mers info
		for(size_t kmer_index=0; kmer_index<m_kmer_temp.size(); ++kmer_index) {
			it_hash = m_kmers_to_index.find(m_kmer_temp[kmer_index]); // Find kmer index
			if(it_hash != m_kmers_to_index.end())  //if in the list
				m_container[(it_hash->second)+hashmap_i] |= or_val; //update container
		}
	}
}

