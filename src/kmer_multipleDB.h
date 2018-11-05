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
#ifndef KMER_MULTIPLEDB
#define KMER_MULTIPLEDB

#include <queue>

#include "kmer_general.h"
#include "kmer_DB.h"

class kmer_heap;
struct bedbim_handle; 

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
		bool load_kmers(const uint64 &batch_size, const kmer_set &set_kmers_to_use);
		inline bool load_kmers(const uint64 &batch_size)
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

		void add_kmers_to_heap(kmer_heap &kmers_and_scores, const std::vector<uint64> &scores, 
				const std::size_t &min_cnt) const;


		inline const std::vector<std::string> get_dbs_names() 
		{ return m_db_names_table; }	// return the list of db names in this class

		// return the indices of DB names inserted in the class DBs


		void clear() {m_kmers.resize(0); m_kmers_table.resize(0);} // clear hashtable

	private:                                
		std::vector<std::string> m_db_names_db_file; 	// Names of DBs in the table file
		std::vector<std::string> m_db_names_table;	// Names of DB to use (table will be squeezed)
		std::size_t m_accessions_db_file;		// Number of DBs in table file 
		std::size_t m_accessions;			// Number of DBs to use (real number of "columns" in save table)
		
		/* tables (file&memory) related parameteres */
		std::ifstream m_kmer_table_file;		// Handel to the table file
		std::size_t m_hash_words_db_file;		// Number of words for "row" in the table file
		std::size_t m_hash_words;			// Number of words for "row" in memory table
		std::vector<uint64> m_kmers;			// List of k-mers (row names)
		std::vector<uint64> m_kmers_table;		// Memort table content
		uint32 m_kmer_len;				// k-mer length
		std::size_t m_left_in_file;			// Number of bytes in table file not read
		std::size_t m_kmer_number;			// Number of k-mers in table file not read
		std::size_t m_kmer_loaded;			// Number of k-mers loaded

		/* Variables for mapping between table in file to saved table in class */
		std::vector<std::size_t> m_map_word_index;	// Word index in table file
		std::vector<std::size_t> m_map_bit_index;	// Bit index in the word

		/* General variables of class */
		bool m_verbose; 				// if we want to output info a log the way

		void create_map_from_all_DBs();
		double calculate_kmer_score(
				const std::size_t kmer_index, 
				const std::vector<uint64> &scores, 
				const double score_sum,
				const uint64 min_in_group = 5) const; 
};

/***********************************************************************************************************/
typedef std::pair<uint64, double> kmer_score;

struct cmp_second
{
	inline bool operator() (const kmer_score& left, const kmer_score& right) const
	{return (left.second*left.second) > (right.second*right.second);} 
};

typedef std::priority_queue<kmer_score, std::vector<kmer_score>, cmp_second> kmer_score_priority_queue;
/**
 * @class kmer_heap
 * @brief save a priority score of kmers
 */
class kmer_heap {
	public:
		kmer_heap(std::size_t max_results);
		void add_kmer(const uint64 &k, const double &score);

		void output_to_file(const std::string &filename) const;
		void output_to_file_with_scores(const std::string &filename) const;
		// it would be nice to some how create a histogram of all the scores along the way...
		// (not only the ones we keep)
		void plot_stat() const;
		inline void empty_heap() {m_best_kmers = kmer_score_priority_queue();} // empty heap content
		kmer_set get_kmer_set() const;	
	private:
		std::size_t m_n_res;
		kmer_score_priority_queue m_best_kmers; // heap that will contain the scores
		//		std::priority_queue<kmer_score, std::vector<kmer_score>, cmp_second> m_best_kmers;
		size_t cnt_kmers;
		size_t cnt_pops;
		size_t cnt_push;	
};



/**
 * @class kmer_multipleDB_merger
 * @brief Taking many kmer_DB's and merge to one table on the disk.
 * Filtering of k-mers can be done in this stage (minimum/maximum count, or part of predefined set)
 */
class kmer_multipleDB_merger {
	public:
		kmer_multipleDB_merger(const std::vector<std::string> &DB_paths,
				const std::vector<std::string> &db_names, 
				const std::string &sorted_kmer_fn,
				const uint32 &kmer_len); // Ctor takes information of the DBs
		
		kmer_multipleDB_merger() = delete; // no default Ctor
		kmer_multipleDB_merger(const kmer_multipleDB_merger&) = delete; // no copy Ctor
		kmer_multipleDB_merger& operator=(const kmer_multipleDB_merger&) = delete; // no equal opertator

		~kmer_multipleDB_merger() {} // desctructor (Dtor of kmer_DB will close the open files)

		// load k-mers from sorted files part of the kmer set
		void load_kmers(const uint64 &iter, const uint64 &total_iter, const kmer_set &set_kmers_to_use);
		inline void load_kmers(const uint64 &iter, const uint64 &total_iter)
		{load_kmers(iter, total_iter, kmer_set());}
		inline void load_kmers() {load_kmers(1ull,1ull);} // load all k-mers in file
		inline void load_kmers(const kmer_set &set_kmers_to_use)
		{load_kmers(1ull, 1ull, set_kmers_to_use);}// load all k-mers in file part of the input kmer set

		void output_to_table(std::ofstream& T) const;
		void output_table_header(std::ofstream& T) const;

		void clear_content(); // clear container 

	private:                                
		std::vector<kmer_DB> m_DBs; // handles to the k-mer files
		std::vector<std::string> m_db_names; // names of the used DBs (accession indices)
		std::vector<uint64> m_kmer_temp; // temp vector to hold read k-mers
		std::size_t m_accessions; // Number of accessions

		std::size_t m_hash_words; // How many words we need to hold the k-mers
	
		my_hash kmers_to_index;
		std::vector<uint64> container;

		uint32 m_kmer_len;

		kmer_DB_sorted_file m_possible_kmers;
};

/**
 * @struct bedbim_handle
 * @brief  holds handles to a bed & bim files (and write the header of the bed file)
 */
struct bedbim_handle {
	bedbim_handle(const std::string &base_name):
		f_bed(base_name + ".bed",std::ios::binary),
		f_bim(base_name + ".bim", std::ios::out) 
	{ f_bed << (char)0x6C << (char)(0x1B) << (char)(0x01);} // header of bed file
	~bedbim_handle() {close();}
	bedbim_handle(bedbim_handle&& o) = default;
	void close();

	std::ofstream f_bed;
	std::ofstream f_bim;
};

#endif

