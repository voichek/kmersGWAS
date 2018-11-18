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

#ifndef KMER_MULTIPLEDB_H
#define KMER_MULTIPLEDB_H

#include "kmer_general.h"
#include "kmers_single_database.h"
#include "best_associations_heap.h"
/**
 * @class kmer_multipleDB
 * @brief class holds the information of presence/absence of k-mers
 * the class should have the ability to contain kmers presence absence information from many  DBs
 * other then that it shouldn't contain other information.
 * Different load options can be available. 
 */
class kmer_multipleDB {
	public:
		kmer_multipleDB(const std::string &merge_db_file, 
				const std::vector<std::string> &db_names,
			        const std::vector<std::string> &db_to_use,	
				const uint32 &kmer_len); // Ctor takes information of the DBs
		kmer_multipleDB() = delete; // no default Ctor
		kmer_multipleDB(const kmer_multipleDB&) = delete; // no copy Ctor
		kmer_multipleDB& operator=(const kmer_multipleDB&) = delete; // no equal opertator

		~kmer_multipleDB() {} // desctructor (Dtor of kmer_DB will close the open files)

		// load k-mers from sorted files part of the kmer set
		bool load_kmers(const uint64_t &batch_size, const kmer_set &set_kmers_to_use);
		inline bool load_kmers(const uint64_t &batch_size)
		{return load_kmers(batch_size, kmer_set());}
		inline bool load_kmers() {return load_kmers(NULL_KEY);} // load all k-mers in file
		inline bool load_kmers(const kmer_set &set_kmers_to_use)
		{return load_kmers(NULL_KEY, set_kmers_to_use);}// load all k-mers in file part of the input kmer set

		void output_kmers_textual() const; // output all k-mers found in the hash to stdout

		// will output a plink bed file along a binary with kmers in the same order
		void output_plink_bed_file(const std::string &base_name) const;
		void output_plink_bed_file(bedbim_handle &f, const kmer_set &set_kmers) const;
		inline void output_plink_bed_file(bedbim_handle &f) const 
		{output_plink_bed_file(f, kmer_set());}
		size_t output_plink_bed_file(bedbim_handle &f, const std::vector<kmer_to_print> &kmer_list, size_t index) const;

		void add_kmers_to_heap(kmer_heap &kmers_and_scores, std::vector<float> scores, 
				const std::size_t &min_cnt) const;


		inline const std::vector<std::string> get_dbs_names() 
		{ return m_db_names_table; }	// return the list of db names in this class

		// return the indices of DB names inserted in the class DBs


		void clear() {m_kmers.resize(0); m_kmers_table.resize(0); m_kmers_popcnt.resize(0);} // clear hashtable

	private:                                
		std::vector<std::string> m_db_names_db_file; 	// Names of DBs in the table file
		std::vector<std::string> m_db_names_table;	// Names of DB to use (table will be squeezed)
		std::size_t m_accessions_db_file;		// Number of DBs in table file 
		std::size_t m_accessions;			// Number of DBs to use (real number of "columns" in save table)
		
		/* tables (file&memory) related parameteres */
		std::ifstream m_kmer_table_file;		// Handel to the table file
		std::size_t m_hash_words_db_file;		// Number of words for "row" in the table file
		std::size_t m_hash_words;			// Number of words for "row" in memory table
		std::vector<uint64_t> m_kmers;			// List of k-mers (row names)
		std::vector<uint64_t> m_kmers_table;		// Memort table content
		std::vector<double> m_kmers_popcnt;		// Memort table content
		uint32 m_kmer_len;				// k-mer length
		std::size_t m_left_in_file;			// Number of bytes in table file not read
		std::size_t m_kmer_number;			// Number of k-mers in table file not read
		std::size_t m_kmer_loaded;			// Number of k-mers loaded
		std::size_t m_row_offset;			// Number of k-mers loaded, not including the current batch
		/* Variables for mapping between table in file to saved table in class */
		std::vector<std::size_t> m_map_word_index;	// Word index in table file
		std::vector<std::size_t> m_map_bit_index;	// Bit index in the word

		/* General variables of class */
		bool m_verbose; 				// if we want to output info a log the way

		void create_map_from_all_DBs();
		double calculate_kmer_score(
				const std::size_t kmer_index, 
				const std::vector<float> &scores, 
				const float score_sum,
				const uint64_t min_in_group = 5) const; 
		void write_PA(const std::string &name, const size_t &kmer_i, bedbim_handle &f) const;
};



#endif

