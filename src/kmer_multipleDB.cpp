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
	m_kmers_pa(5000000),
	m_hash_words((m_accessions+64-1)/64)
//	m_kmers_pa((m_accessions+64-1)/64) // a way to divide and take ceil
{
	// build all the kmer_DB objects and open the sorted k-mer file
	for(size_t i=0; i < db_names.size(); i++) {
		m_DBs.emplace_back(path_to_DBs +"/" + db_names[i], db_names[i]);
		m_DBs.back().open_sorted_kmer_file(sorted_kmer_fn);
	}
	cerr << "We have " << m_accessions << " accessions and " << m_kmers_pa.size() << " hash-map\n";

	// Initializing the hash maps (for every 64 accessions we have 1)
//	for(size_t i=0; i < m_kmers_pa.size(); i++)
//		m_kmers_pa[i].set_empty_key(-1);
	m_kmers_pa.set_empty_key(0xFFFFFFFFFFFFFFFF);
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
//	for(size_t i=0; i < m_kmers_pa.size(); i++)
//		m_kmers_pa[i].clear();
	m_kmers_pa.clear();

	// To make sure all my hash maps are aligned I will initialize to 0 all the other ones
	vector<uint64> new_kmers;	
	my_multi_hash::iterator it_hash;
	for(size_t acc_i = 0; acc_i < m_accessions; ++acc_i) {
		//cerr << "loading " << acc_i ;
		size_t hashmap_i = acc_i / 64;
		size_t bit_i = acc_i % 64;

		uint64 or_val = 1ull << bit_i;
		array<uint64, WORD64HASHT> new_bits;
		new_bits[0] = 0ull;
		new_bits[1] = 0ull;
		new_bits[hashmap_i] = or_val;
//		cerr << "\tnew bit: "<<bitset<64>(or_val) << endl;
		
		// Reading new-kmers from file
		// new_kmers.resize(0); - the inner function is allready emptying this vector
		m_DBs[acc_i].read_sorted_kmers(m_kmer_temp, current_threshold);
		for(vector<uint64>::iterator it = m_kmer_temp.begin(); it != m_kmer_temp.end(); ++it) {
//			it_hash = m_kmers_pa[hashmap_i].find(*it);
			it_hash = m_kmers_pa.find(*it);			
			if(it_hash == m_kmers_pa.end()) 
			{ 

				m_kmers_pa.insert(my_multi_hash::value_type(*it, new_bits));
	//			new_kmers.push_back(*it);
			}
			else {
				it_hash->second[hashmap_i] |= or_val;
			}
		}
		// Add the new k-mers to the other hash_maps
	}
}

void kmer_multipleDB::plot_textual_hash_map(const std::vector<double> &phenotypes) {
//	my_multi_hash::iterator it_hash;
//	for(auto it : m_kmers_pa) {
	for(my_multi_hash::iterator it=m_kmers_pa.begin(); it != m_kmers_pa.end(); ++it) {		
		double x = calculate_kmer_score(it, phenotypes);
		if (x*x > 16.) { 
			cout << bitset<64>(it->first); // << "\t" << bitset<64>(it.second);
			for(size_t i=0; i<m_hash_words; i++) 
			{
				//			it_hash = m_kmers_pa[i].find(it.first);
				cout << "\t" << bitset<64>(it->second[i]);
			}
			cout << "\t" << x << endl;
		}
		//		cout << x << endl;
	}
}

double kmer_multipleDB::calculate_kmer_score(my_multi_hash::iterator& it, 
		const std::vector<double> &scores,
		double min_in_group) {
	double Ex0, Ex1, E2x0, E2x1, Sx0, Sx1, N0, N1, bit;
	Ex0 = Ex1 = Sx0 = Sx1 = N0 = N1 = E2x0 = E2x1 = 0;
	for(size_t i=0; i<m_accessions; i++) {
		size_t hashmap_i = i / 64;
		size_t bit_i = i % 64;
		bit = double((it->second[hashmap_i] >> bit_i)&1);
		Ex1 += bit*scores[i];
		E2x1 += bit*scores[i]*scores[i];
		Ex0 += (1-bit)*scores[i];
		E2x0 += (1-bit)*scores[i]*scores[i];
		N1 = N1+bit;
	}
	N0 = m_accessions - N1;
	//	cout << "\t" << N0 << "\t" << N1;
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

		//		students_t dist(v);
		//		double q = cdf(complement(dist, fabs(t_stat)));
		//		return t_test;
	}
	else {return 0;}

	// Degrees of freedom:
}
