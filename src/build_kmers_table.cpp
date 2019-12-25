///
///      @file  build_kmers_table.cpp
///     @brief  This program will create a table of presence absence accross DBs
//		
//		Taking the sorted general kmers list as well as the sorted kmers from 
//		each accessions, create a table with kmers as rows and accessions as collumns
//		where each cell indicate if the kmers was found in the accessions or not
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  10/29/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include <cxxopts/include/cxxopts.hpp>

#include "kmer_general.h"
#include "kmers_merge_multiple_databaes.h"

using namespace std;

int main(int argc, char *argv[]) {
	cxxopts::Options options("build_kmers_table", "Build the k-mers table");
	try
	{
		options.add_options()
			("l,list_kmers_files", "list of seperate k-mers files", cxxopts::value<string>())
			("k,kmers_len", "length of k-mers", cxxopts::value<size_t>())
			("a,all_kmers", "path to file with all k-mers", cxxopts::value<string>())
			("o,output", "prefix for kmers-table files", cxxopts::value<string>())
			("help", "print help")
			;
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cerr << options.help() << endl;
			exit(0);
		}
		vector<string> required_parametrs({"list_kmers_files","kmers_len","all_kmers","output"});
		for(size_t i=0; i<required_parametrs.size(); i++) {
			if(result.count(required_parametrs[i]) == 0) {
				cerr << required_parametrs[i] << " is a required parameter" << endl;
				cerr << options.help() << endl;
				exit(1);
			}
		}
		// Load parameters
		string fn_list_kmers_files(result["list_kmers_files"].as<string>());
		string fn_all_kmers(result["all_kmers"].as<string>());
		size_t kmer_len(result["kmers_len"].as<size_t>());
		string output_base(result["output"].as<string>());

		// Check if all input files exist
		vector<string> required_files({fn_list_kmers_files,fn_all_kmers});
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

		// Read accessions to use
		vector<AccessionPath> kmers_handles = read_accessions_path_list(fn_list_kmers_files);
		vector<string> kmers_filenames, accessions_names;

		ofstream fout_names(output_base + ".names", ios::binary);
		for(size_t i=0; i<kmers_handles.size(); i++) {
			accessions_names.push_back(kmers_handles[i].name);
			fout_names << kmers_handles[i].name << endl;

			if(!is_file_exist(kmers_handles[i].path)) {
				cerr << "Couldn't find file: " << kmers_handles[i].path << endl;
				exit(1);
			}
			kmers_filenames.push_back(kmers_handles[i].path);
		}
		fout_names.close();
		cerr << "Create merger" << endl;	
		MultipleKmersDataBasesMerger merger(kmers_filenames, accessions_names, fn_all_kmers, kmer_len);

		cerr << "Opens file" << endl;
		ofstream fout(output_base + ".table", ios::binary);
		merger.output_table_header(fout);
		uint64_t total_iter = 5000; // Can move this to the input parameters if needed
		for(uint i=1; i<=(total_iter+1); i++) { // +1 is for debugging - should be empty
			cerr << i << " / " << total_iter << " : Loading k-mers" << endl;
			merger.load_kmers(i, total_iter);
			merger.output_to_table(fout);
		}	
		cerr << "close file" << endl;
		fout.close();
	}  catch (const cxxopts::OptionException& e)
	{
		cerr << "error parsing options: " << e.what() << endl;
		cerr << options.help() << endl;
		exit(1);
	}
	return 0;
}
