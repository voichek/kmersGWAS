#ifndef KMER_DB_H
#define KMER_DB_H

#include "kmer_general.h"


/**
 * @class kmer_DB_sorted_file
 * @brief managed accesses to a sorted k-mer file
 */

class kmer_DB_sorted_file {
	public:
	kmer_DB_sorted_file(); // Default constructor
	kmer_DB_sorted_file(const std::string &filename);
	kmer_DB_sorted_file(const kmer_DB_sorted_file&) = delete;
	kmer_DB_sorted_file(kmer_DB_sorted_file&&) = default;
	~kmer_DB_sorted_file(); // Make sure the file is close

	void open_file(const std::string &filename);
	void close_file();
	
	void load_kmers_upto_x(const uint64 &threshold, std::vector<uint64> &kmers);

	uint64 get_kmer_count(); // get the number of kmers in the file

	private:
	void init();
	void read_kmer();

	// variabled for reading sortet kmer file
	ifstream m_fin;			// handle to file
	uint64 m_last_kmer;		// last kmer loaded
	uint64 m_kmers_in_file; // number of k-mers in file (calcluated once in init)
	uint64 m_kmers_count;			// last kmer index in file
};



/* class Kmer_DB will represent a database of kmers
 * First, it will have a handel on the DB created by KMC and can easily open and iterate over it
 * Second, it will be able to intersect the DB with a set of kmers and output the result to an output file
 */
class Kmer_DB {
	public:
	Kmer_DB(const std::string dir_path, const std::string db_name); // Constructor
	Kmer_DB(const Kmer_DB &x) = delete;
	Kmer_DB(Kmer_DB&&) = default;
	~Kmer_DB() {}; //Destructor (compiler generate one automatically, but safe to be explicit)

	// return the name of the DB
	std::string get_name() { return m_db_name; }
	
	// return the directory in which the KMC DB is found
	std::string get_dir_path() { return m_dir_path; }

	// go over the KMC DB and output to file only the kmers in the kmers set inputed
	// output to file should be numerically ordered
	void intersect_kmers(const kmer_set& kmers_to_use, std::string file_name);

	// open kmer sorted file name for reading
	void open_sorted_kmer_file(const std::string& filename);

	// reads all k-mers until getting to somethreshold
	void read_sorted_kmers(std::vector<uint64> &kmers, uint64 threshold = 0xFFFFFFFFFFFFFFFF);

	// Counts how many times each k-mer appeared and plot to std::cout
	std::vector<std::size_t> calculate_kmers_counts_histogram();

	// get a hanle to the KMC DB to run over its k-mers 
	CKMCFile get_KMC_handle();

	private:
//	Kmer_DB(const Kmer_DB&); // Only declare and no implementation - to avoid default by compiler
	Kmer_DB& operator=(const Kmer_DB&) = delete; //same for the = operator


	std::string m_db_name;
	std::string m_dir_path;
	
	kmer_DB_sorted_file m_sorted_kmers_f;
};



#endif

