/**
 *       @file  kmers_merge_multiple_databaes.h
 *      @brief  
 *
 *	The merge kmers table is NOT promised to have k-mers ordered!
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

#ifndef KMERS_MERGE_MULTIPLE_DATABASES 
#define KMERS_MERGE_MULTIPLE_DATABASES 

#include "kmer_general.h"
#include "kmers_single_database.h"
/**
 * @class MultipleKmersDataBasesMerger
 * @brief Taking many KmersSingleDataBase's and merge to one table on the disk.
 * Filtering of k-mers can be done in this stage (minimum/maximum count, or part of predefined set)
 */
class MultipleKmersDataBasesMerger {
	public:
		MultipleKmersDataBasesMerger(const std::vector<std::string> &DB_paths,
				const std::vector<std::string> &db_names, 
				const std::string &sorted_kmer_fn,
				const uint32 &kmer_len); // Ctor takes information of the DBs
		
		MultipleKmersDataBasesMerger() = delete; // no default Ctor
		MultipleKmersDataBasesMerger(const MultipleKmersDataBasesMerger&) = delete; // no copy Ctor
		MultipleKmersDataBasesMerger& operator=(const MultipleKmersDataBasesMerger&) = delete; // no equal opertator

		~MultipleKmersDataBasesMerger() {} // desctructor (Dtor of KmersSingleDataBase will close the open files)

		// load k-mers from sorted files part of the kmer set
		void load_kmers(const uint64_t &iter, const uint64_t &total_iter, const KmersSet &set_kmers_to_use);
		inline void load_kmers(const uint64_t &iter, const uint64_t &total_iter)
		{load_kmers(iter, total_iter, KmersSet());}
		inline void load_kmers() {load_kmers(1ull,1ull);} // load all k-mers in file
		inline void load_kmers(const KmersSet &set_kmers_to_use)
		{load_kmers(1ull, 1ull, set_kmers_to_use);}// load all k-mers in file part of the input kmer set

		void output_to_table(std::ofstream& T) const;
		void output_table_header(std::ofstream& T) const;

		void clear_content(); // clear container 

	private:                                
		std::vector<KmersSingleDataBase> m_DBs; // handles to the k-mer files
		std::vector<std::string> m_db_names; // names of the used DBs (accession indices)
		std::vector<uint64_t> m_kmer_temp; // temp vector to hold read k-mers
		std::size_t m_accessions; // Number of accessions

		std::size_t m_hash_words; // How many words we need to hold the k-mers
	
		KmerUint64Hash kmers_to_index;
		std::vector<uint64_t> container;

		uint32 m_kmer_len;

		KmersSingleDataBaseSortedFile m_possible_kmers;
};


#endif

