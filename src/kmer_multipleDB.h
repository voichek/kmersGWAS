///
///      @file  kmer_multipleDB.h
///     @brief  Decleration of kmer_multiDB class
///
/// Class explained bellow will combine the k-mers from multiple k-mer DBs so tests
/// on comparing the presence of k-mers in a few DBs can be compared
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
#include <array>

#include "kmer_general.h"
#include "kmer_DB.h"

#define WORD64HASHT 16

typedef google::dense_hash_map<uint64, std::array<uint64, WORD64HASHT>, Hash64> my_multi_hash; 
/**
 * @class kmer_multipleDB
 * @brief 
 */

class kmer_multipleDB {
	public:
		kmer_multipleDB(const string& path_to_DBs, 
				const std::vector<std::string>& db_names, 
				const string& sorted_kmer_fn);
		~kmer_multipleDB() {};
		kmer_multipleDB() = delete;
		kmer_multipleDB(const kmer_multipleDB&) = delete;
		kmer_multipleDB& operator=(const kmer_multipleDB&) = delete;
		
		void load_kmers(const uint64 &iter, const uint64 &total_iter);
		void plot_textual_hash_map(const std::vector<double> &phenotypes);
		double calculate_kmer_score(my_multi_hash::iterator& it, const std::vector<double> &scores,
		double min_in_group = 5.);
	private:
		std::vector<Kmer_DB> m_DBs;
		std::vector<uint64> m_kmer_temp;
		std::size_t m_accessions;
//		std::vector<my_hash> m_kmers_pa;		
		my_multi_hash m_kmers_pa;		
		size_t m_hash_words;
};
