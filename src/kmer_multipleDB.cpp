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
#include <bitset>
using namespace std;

///
/// @brief  Ctor of kmer_multipleDB - 
/// @param  path_to_DBs:	a path in the filesystem to directories containing the files ok k-mers
//			db_names:		a list of names of all the DBs to use
//			sorted_kmer_fn:	filename inside each subdirectory containing the sorted k-mer list
/// @return 
///
kmer_multipleDB::kmer_multipleDB(const string& path_to_DBs, const std::vector<std::string>& db_names, 
		const string& sorted_kmer_fn):
	m_DBs(),
	m_kmer_temp(),
	m_accessions(db_names.size()),
	m_kmers_pa((m_accessions+64-1)/64) // a way to divide and take ceil
{
	// build all the kmer_DB objects and open the sorted k-mer file
	for(size_t i=0; i < db_names.size(); i++) {
		m_DBs.emplace_back(path_to_DBs +"/" + db_names[i], db_names[i]);
		m_DBs.back().open_sorted_kmer_file(sorted_kmer_fn);
	}
	cerr << "We have " << m_accessions << " accessions and " << m_kmers_pa.size() << " hash-map\n";

	// Initializing the hash maps (for every 64 accessions we have 1)
	for(size_t i=0; i < m_kmers_pa.size(); i++)
		m_kmers_pa[i].set_empty_key(-1);

}


/**
 * @brief   kmer_multipleDB::load_kmers - load a subset of k-mers from sorted files
 * @param   Which iteration (iter) is it from all iterations (total_iter)
 * @return  
 */
void kmer_multipleDB::load_kmers(const uint64 &iter, const uint64 &total_iter) {
	// as k-mers are 31 bp - the largest possible value is 0011111111...1111
	uint64 current_threshold = ((0x3FFFFFFFFFFFFFFF / total_iter)+1)*iter;
	cerr << "threshold: " << bitset<64>(current_threshold) << endl; 
	// Clear the content of the hash maps
	for(size_t i=0; i < m_kmers_pa.size(); i++)
		m_kmers_pa[i].clear();

	// To make sure all my hash maps are aligned I will initialize to 0 all the other ones
	vector<uint64> new_kmers;	
	my_hash::iterator it_hash;
	for(size_t acc_i = 0; acc_i < m_accessions; ++acc_i) {
		cerr << "loading " << acc_i ;
		size_t hashmap_i = acc_i / 64;
		size_t bit_i = acc_i % 64;

		uint64 or_val = 1ull << bit_i;
		cerr << "\tnew bit: "<<bitset<64>(or_val) << endl;
		
		// Reading new-kmers from file
		// new_kmers.resize(0); - the inner function is allready emptying this vector
		m_DBs[acc_i].read_sorted_kmers(m_kmer_temp, current_threshold);
		for(vector<uint64>::iterator it = m_kmer_temp.begin(); it != m_kmer_temp.end(); ++it) {
			it_hash = m_kmers_pa[hashmap_i].find(*it);
			if(it_hash == m_kmers_pa[hashmap_i].end()) 
			{ 
				m_kmers_pa[hashmap_i].insert(my_hash::value_type(*it, or_val));
	//			new_kmers.push_back(*it);
			}
			else {
				it_hash->second |= or_val;
			}
		}
		// Add the new k-mers to the other hash_maps
	//	for(size_t i=0; i < m_kmers_pa.size(); i++) {
	//		if(i != hashmap_i) {
	//			cerr << "adding new to:" << i << " : " << new_kmers.size() << endl;
	//			for(vector<uint64>::iterator it = new_kmers.begin(); it != new_kmers.end(); ++it) {
	//				m_kmers_pa[i].insert(my_hash::value_type(*it,0ull));
	//			}
	//		}
	//	}
	}
}

void kmer_multipleDB::plot_textual_hash_map() {
	my_hash::iterator it_hash;
	for(auto it : m_kmers_pa[0]) {
		cout << bitset<64>(it.first) << "\t" << bitset<64>(it.second);
		for(size_t i=1; i<m_kmers_pa.size(); i++) 
		{
			it_hash = m_kmers_pa[i].find(it.first);
			cout << "\t" << bitset<64>(it_hash->second);
		}
		cout << endl;
	}
}
