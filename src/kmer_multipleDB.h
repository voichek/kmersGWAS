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

#include <array>
#include <queue>

#include "kmer_general.h"
#include "kmer_DB.h"

#define WORD64HASHT 16

typedef google::dense_hash_map<uint64, std::array<uint64, WORD64HASHT>, Hash64> my_multi_hash; 
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
		kmer_multipleDB(const std::vector<string> &DB_paths, // Intentionally wrote string and not std:: as the KMC package include namespace std and I want to fix it 
				const std::vector<std::string> &db_names, 
				const std::string &sorted_kmer_fn); // Ctor takes information of the DBs
		kmer_multipleDB() = delete; // no default Ctor
		kmer_multipleDB(const kmer_multipleDB&) = delete; // no copy Ctor

		kmer_multipleDB& operator=(const kmer_multipleDB&) = delete; // no equal opertator

		~kmer_multipleDB() {} // desctructor (Dtor of kmer_DB will close the open files)

		// load k-mers from sorted files part of the kmer set
		void load_kmers(const uint64 &iter, const uint64 &total_iter, const kmer_set &set_kmers_to_use);
		inline void load_kmers(const uint64 &iter, const uint64 &total_iter)
		{load_kmers(iter, total_iter, kmer_set());}
		inline void load_kmers() {load_kmers(1ull,1ull);} // load all k-mers in file
		inline void load_kmers(const kmer_set &set_kmers_to_use)
		{load_kmers(1ull, 1ull, set_kmers_to_use);}// load all k-mers in file part of the input kmer set

		void output_kmers_textual() const; // output all k-mers found in the hash to stdout
		void output_kmers_binary(const std::string &filename) const;

		// will output a plink bed file along a binary with kmers in the same order
		void output_plink_bed_file(const std::string &base_name) const;
		void output_plink_bed_file(bedbim_handle &f, const kmer_set &set_kmers) const;
		inline void output_plink_bed_file(bedbim_handle &f) const 
		{output_plink_bed_file(f, kmer_set());}

		// Continuous phenotype (e.g flowering time 1,2,3..100)
		void add_kmers_to_heap(kmer_heap &kmers_and_scores, const std::vector<double> &scores, 
				const std::vector<std::string> &names_scores) const;
		void add_kmers_to_heap(kmer_heap &kmers_and_scores, const std::vector<double> &scores) const; //efficent procedure
		// Categorical phenotype (e.g resistence yes/no) - never debugged!
		void add_kmers_to_heap(kmer_heap &kmers_and_scores, const std::vector<uint64> &scores, 
				const std::size_t &min_cnt) const;


		inline const std::vector<std::string> get_dbs_names() 
		{ return m_db_names; }	// return the list of db names in this class

		// return the indices of DB names inserted in the class DBs
		std::vector<std::size_t> get_dbs_indices(const std::vector<std::string> &names) const;

		inline std::size_t get_hashtable_size() const {return m_kmers_pa.size();}

		void clear_hashtable() {m_kmers_pa.clear();} // clear hashtable

	private:                                
		std::vector<kmer_DB> m_DBs; // handles to the k-mer files
		std::vector<std::string> m_db_names; // names of the used DBs (accession indices)
		std::vector<uint64> m_kmer_temp; // temp vector to hold read k-mers
		std::size_t m_accessions; // Number of accessions
		my_multi_hash m_kmers_pa; // hash table where k-mers are combines
		std::size_t m_hash_words; // How many words we need to hold the k-mers (though we use a const amount)
		bool m_verbose; // if we want to output info a log the way

		// Functions that calculate association scores
		// 1. Two sided t-test
		double calculate_kmer_score(const my_multi_hash::const_iterator& it, 
				const std::vector<double> &scores, const std::vector<std::size_t> &accession_index,
				const double min_in_group = 5.) const;
		double calculate_kmer_score(const my_multi_hash::const_iterator& it, 
				const vector<double> &scores, 				const vector<double> &scores2,
				const double score_sum,			const double score2_sum,
				const double min_in_group =5.) const;
		double calculate_kmer_score(const my_multi_hash::const_iterator& it, 
				const vector<uint64> &scores,	const vector<uint64> &scores2,
				const uint64 score_sum,			const uint64 score2_sum,
				const uint64 min_in_group =5) const;

		// 2. fisher exact test if we have only two options for scores (e.g resistence +/-)
		double calculate_kmer_score(const my_multi_hash::const_iterator& it, 
				const std::vector<size_t> &scores, const std::vector<std::size_t> &accession_index,
				double min_in_group = 5.) const;

		// 3. Wilcoxon-rank-sum-test - in case we want to be consistent with Nordborg 107 phenotype Nature
		// double kmer_score_wilcoxon_rank_sum_test(const my_multi_hash::iterator& it, const std::vector<double> &scores,
		// double min_in_group = 5.)
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
 * @struct bedbim_handle
 * @brief  holds handles to a bed & bim files (and write the header of the bed file)
 */
struct bedbim_handle {
	bedbim_handle(const std::string &base_name):
		f_bed(base_name + ".bed",ios::binary),
		f_bim(base_name + ".bim", ios::out) 
	{ f_bed << (char)0x6C << (char)(0x1B) << (char)(0x01);} // header of bed file
	~bedbim_handle() {close();}
	bedbim_handle(bedbim_handle&& o) = default;
	void close();

	ofstream f_bed;
	ofstream f_bim;
};

#endif

