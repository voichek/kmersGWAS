///
///      @file  histogram_KMC_kmers_counts.cpp 
///     @brief  Create an histogram of k-mers count for KMC DB 
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  08/14/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include <cxxopts/include/cxxopts.hpp>
#include "kmers_single_database.h"

using namespace std;
int main(int argc, char *argv[]) {
	cxxopts::Options options("histogram_KMC_kmers_counts", "Calculate the histogram of k-mers counts from a KMC DB (output to stdout)");
	try
	{
		options.add_options()
			("d,dir", "path to directory with KMC DBs", cxxopts::value<string>())
			("n,name", "prefix name of KMC DB files", cxxopts::value<string>())
			("k,kmers_len", "length of k-mers", cxxopts::value<size_t>())
			("help", "print help")
			;
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cerr << options.help() << endl;
			exit(0);
		}
		vector<string> required_parametrs({"dir", "kmers_len", "name"});
		for(size_t i=0; i<required_parametrs.size(); i++) {
			if(result.count(required_parametrs[i]) == 0) {
				cerr << required_parametrs[i] << " is a required parameter" << endl;
				cerr << options.help() << endl;
				exit(1);
			}
		}
		// Load parameters
		string dir_path(result["dir"].as<string>());
		string db_name(result["name"].as<string>());
		size_t kmer_len(result["kmers_len"].as<size_t>());

		// Check if all input files exist
		vector<string> required_files({dir_path + "/" + db_name + ".kmc_pre", dir_path + "/" + db_name + ".kmc_suf"});
		for(size_t i=0; i<required_files.size(); i++) {
			if(!is_file_exist(required_files[i])) {
				cerr << "Couldn't find file: " << required_files[i] << endl;
				exit(1);
			}
		}
		if((kmer_len > 31) || (kmer_len <10)) {
			cerr << "kmer length has to be between 10-31" << endl;
			exit(1);
		}
		/****************************************************************************************************/
		/* END of parsing and checking input parameters */
		/****************************************************************************************************/
		KmersSingleDataBase acc_kmers(dir_path, db_name, kmer_len);
		vector<size_t> counters = acc_kmers.calculate_kmers_counts_histogram();

		cout << "appearance\tcount" << endl;
		for(size_t i=0; i<counters.size(); i++) 
			cout << i << "\t" << counters[i] << endl;											
	} catch (const cxxopts::OptionException& e)
	{
		cerr << "error parsing options: " << e.what() << endl;
		cerr << options.help() << endl;
		exit(1);
	}
	return 0;
}

