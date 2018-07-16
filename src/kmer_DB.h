#ifndef KMER_DB_H
#define KMER_DB_H

#include "kmer_general.h"
#include <string>

/* class Kmer_DB will represent a database of kmers
 * First, it will have a handel on the DB created by KMC and can easily open and iterate over it
 * Second, it will be able to intersect the DB with a set of kmers and output the result to an output file
 */
class Kmer_DB {
	public:
	Kmer_DB(std::string dir_path, std::string db_name); // Constructor
	
	// return the name of the DB
	std::string get_name() { return m_db_name; }
	
	// return the directory in which the KMC DB is found
	std::string get_dir_path() { return m_dir_path; }

	// go over the KMC DB and output to file only the kmers in the kmers set inputed
	// output to file should be numerically ordered
	void intersect_kmers(const kmer_set& kmers_to_use, std::string file_name);

	// Counts how many times each k-mer appeared and plot to std::cout
	void calculate_kmers_counts_histogram();

	// get a hanle to the KMC DB to run over its k-mers 
	CKMCFile get_KMC_handle();

	private:
	std::string m_db_name;
	std::string m_dir_path;
};

#endif

