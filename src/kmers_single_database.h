///
///      @file  kmers_single_database.h
///     @brief  Decleration of class KmersSingleDataBase
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  12/01/19
///   Revision  $Id: doxygen.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2019, Yoav Voichek
///
///This source code is released for free distribution under the terms of the
///GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///


#ifndef KMER_DB_H
#define KMER_DB_H

#include "kmer_general.h"


/**
 * @class KmersSingleDataBaseSortedFile
 * @brief Managed accesses to a sorted k-mer file
 */

class KmersSingleDataBaseSortedFile {
	public:
	KmersSingleDataBaseSortedFile(); // Default constructor
	KmersSingleDataBaseSortedFile(const std::string &filename);
	KmersSingleDataBaseSortedFile(const KmersSingleDataBaseSortedFile&) = delete;
	KmersSingleDataBaseSortedFile(KmersSingleDataBaseSortedFile&&) = default;
	~KmersSingleDataBaseSortedFile(); // Make sure the file is close

	void open_file(const std::string &filename);
	void close_file();
	
	void load_kmers_upto_x(const uint64_t &threshold, std::vector<uint64_t> &kmers);
	void load_kmers_upto_x(const uint64_t &threshold, std::vector<uint64_t> &kmers, std::vector<uint64_t> &flags); 

	uint64_t get_kmer_count() {return m_kmers_in_file;} // get the number of kmers in the file

	private:
	void init();
	void read_kmer();

	// variabled for reading sorted kmer file
	std::ifstream m_fin;		// handle to file
	uint64_t m_last_kmer;		// last kmer loaded
	uint64_t m_flag;		// last flag (2 MSB in kmer)
	uint64_t m_kmers_in_file; 	// number of k-mers in file (calcluated once in init)
	uint64_t m_kmers_count;		// last kmer index in file
};



/**
 * @class KmersSingleDataBase 
 * @bried Represent a database of kmers
 * 1. Handles a DB created by KMC and can easily open and iterate over it
 * 2. Handles a sorted kmers list from the same folder
 * 3. Can intersect the DB with a set of kmers and output the result to an output file
 */
class KmersSingleDataBase {
	public:
	KmersSingleDataBase(const std::string& dir_path, const std::string& db_name, const uint32& kmer_length); // Constructor
	KmersSingleDataBase(const KmersSingleDataBase &x) = delete;
	KmersSingleDataBase(KmersSingleDataBase&&) = default;
	KmersSingleDataBase& operator=(const KmersSingleDataBase&) = delete; //same for the = operator
	~KmersSingleDataBase() {}; //Destructor (compiler generate one automatically, but safe to be explicit)

	// return the name of the DB
	std::string get_name() { return m_db_name; }
	
	// return the directory in which the KMC DB is found
	std::string get_dir_path() { return m_dir_path; }

	// go over the KMC DB and output to file only the kmers in the kmers set inputed
	// output to file should be numerically ordered
	void intersect_kmers(const KmersSet& kmers_to_use, std::string file_name);

	// open kmer sorted file name for reading
	void open_sorted_kmer_file(const std::string& filename) {
		m_sorted_kmers_f.open_file(m_dir_path + "/" + filename);}

	// reads all k-mers until getting to somethreshold
	void read_sorted_kmers(std::vector<uint64_t> &kmers, uint64_t threshold = 0xFFFFFFFFFFFFFFFF) {
		m_sorted_kmers_f.load_kmers_upto_x(threshold, kmers); }
	// Counts how many times each k-mer appeared and plot to std::cout
	std::vector<std::size_t> calculate_kmers_counts_histogram();

	// get a hanle to the KMC DB to run over its k-mers 
	CKMCFile get_KMC_handle();

	private:
	std::string m_db_name;
	std::string m_dir_path;
	KmersSingleDataBaseSortedFile m_sorted_kmers_f;
	uint32 m_kmer_len; 
};

#endif

