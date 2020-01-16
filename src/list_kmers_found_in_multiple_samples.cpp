///
///      @file  list_kmers_found_in_multiple_samples.cpp
///     @brief  Go over a list of k-mer DBs and output only k-mers appearing in at least 
/// 		N (defined by user) DBs. The filtered k-mers are output sorted in a binary format. 
///
//		Read all the sorted kmers files and collect statistics on them:
//		1. In how many accessions they were found
//		2. In how many they apeard in canonized/non-canonized form or both
//
//		Output the general statistics on how many times they were found in
//		which form.
//
//		Filter kmers to the ones apearing at least N times and found in
//		every form at least in some defined percent.
//
//		For example: a kmer K will pass the threshold if it:
//		1. Found in at least 5 acessions
//		2. Found at least 20% of the time in every form
//		*  So if it apeared in 100 accessions, we should 20 accessions
//		with each form.
//
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///    Created  08/14/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmer_general.h"
#include "kmers_single_database.h"

#include <cmath>

#include <cxxopts/include/cxxopts.hpp>

using namespace std;

void output_uint64_t_matrix(const string& fn, const vector<vector<uint64_t> >& M) {
	ofstream fout(fn);
	for(size_t i=0; i<M.size(); i++) {
		for(size_t j=0; j<M[i].size(); j++) {
			if(j>0)
				fout << "\t";
			fout <<  M[i][j];
		}
		fout << endl; 
	}
	fout.close();
}

int main(int argc, char *argv[]) {
	cxxopts::Options options("list_kmers_found_in_multiple_samples", "Combines and filters information from all samples k-mers lists to one sorted k-mers list");
	try
	{
		options.add_options()
			("l,list_kmers_files", "list of separate k-mers files", cxxopts::value<string>())
			("k,kmers_len", "length of k-mers", cxxopts::value<size_t>())
			("mac", "minor allele count (minimum allowed appearence of a k-mer)", cxxopts::value<size_t>())
			("p,min_strand_percent", "minimum percent of apperence in each strand form", cxxopts::value<double>())
			("o,output", "path to output file", cxxopts::value<string>())
			("help", "print help")
			;
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cerr << options.help() << endl;
			exit(0);
		}
		vector<string> required_parametrs({"list_kmers_files", "mac",  "kmers_len", "min_strand_percent", "output"});
		for(size_t i=0; i<required_parametrs.size(); i++) {
			if(result.count(required_parametrs[i]) == 0) {
				cerr << required_parametrs[i] << " is a required parameter" << endl;
				cerr << options.help() << endl;
				exit(1);
			}
		}

		// Load parameters
		string fn_kmers_list(result["list_kmers_files"].as<string>());
		size_t minimum_kmer_count(result["mac"].as<size_t>());
		size_t kmer_len(result["kmers_len"].as<size_t>());
		double minimum_strand_per(result["min_strand_percent"].as<double>());
		string fn_output(result["output"].as<string>());

		// Check if all input files exist
		vector<string> required_files({fn_kmers_list});
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
		vector<AccessionPath> sorted_kmers_fn = read_accessions_path_list(fn_kmers_list);
		vector<KmersSingleDataBaseSortedFile> list_handles;

		// Open all kmers files
		size_t N = sorted_kmers_fn.size();
		for(size_t i=0; i<N; i++) { 
			if(!is_file_exist(sorted_kmers_fn[i].path)) {
				cerr << "Couldn't find file: " << sorted_kmers_fn[i].path << endl;
				exit(1);
			}
			list_handles.emplace_back(sorted_kmers_fn[i].path);
		}

		// Build the hash table that will contain all the kmers counts
		KmerUint64Hash kmers_hash(10000000); // Initial hash_table_size from user
		kmers_hash.set_empty_key(NULL_KEY);

		// Defining variables to use while going over the k-mers
		vector<uint64_t> kmers, flags, unique_kmers;
		KmerUint64Hash::iterator it_hash;

		// statistics collection
		vector<uint64_t> shareness(N+1, 0);	
		vector<vector<uint64_t> > only_canonical(N+1,vector<uint64_t>(N+1,0));
		vector<vector<uint64_t> > only_non_canonical(N+1,vector<uint64_t>(N+1,0));
		vector<vector<uint64_t> > both_forms(N+1,vector<uint64_t>(N+1,0));

		// We want to keep track for each kmer if it apeared in canonical/non-canonical or both
		// as the counters have 64 bits and we assume we will have less than 1,000,000 individuals
		// we can put 3 counters in the same word. The following array define the constant to add.
		uint64_t adders[] = {1ull+(1ull<<20ull), 1ull+(1ull<<40ull), 1ull+0ull}; 

		// Output files:
		ofstream file_out_kmers(fn_output, ios::binary); // kmers passing filters
		ofstream file_non_pass_kmers(fn_output + ".no_pass_kmers"); // kmers not passing filters
		file_non_pass_kmers << "kmer\tcount_all\tcanonical\tnon-canonical\tboth" << endl;

		uint64_t STEPS = 5000;
		uint64_t cnt_no_pass(0), cnt_pass(0), cnt_low_MAC;
		for(uint64_t step_i=1; step_i<=(STEPS+1); step_i++) { // the +1 in STEPS is for debugging
			kmers_hash.clear(); // Empty the hash table
			unique_kmers.resize(0); // Empty unique kmers collection

			// Set threshold until which kmer to read
			uint64_t current_threshold = ((((1ull << (kmer_len*2ull))-1ull) / STEPS)+1) * step_i;
			cerr << step_i << " / " << STEPS << "\t:\t" << bitset<WLEN>(current_threshold);

			// Go over all the kmers and collect information
			for(size_t i=0; i<N; i++) { //over individuals
				list_handles[i].load_kmers_upto_x(current_threshold, kmers, flags); // read kmers & flags
				for(size_t kmer_i=0; kmer_i < kmers.size(); kmer_i++) {
					uint64_t cur_adder = adders[flags[kmer_i]-1]; // Flag should never be zero!

					it_hash = kmers_hash.find(kmers[kmer_i]);
					if(it_hash == kmers_hash.end()) { // new kmer
						kmers_hash.insert(KmerUint64Hash::value_type(kmers[kmer_i], cur_adder));
						unique_kmers.push_back(kmers[kmer_i]);
					} else
						it_hash->second += cur_adder;
				}
			}
			cerr << "\tunique kmers: " << unique_kmers.size() << endl;
			sort(unique_kmers.begin(), unique_kmers.end()); // Sort so the end list will also be sorted
			// output kmers that pass threshold and collect statistics
			for(size_t i=0; i<unique_kmers.size(); i++) {
				it_hash = kmers_hash.find(unique_kmers[i]); // kmers from unique_kmers have to be in hash
				uint64_t counts = it_hash->second;
				uint64_t count_all       = counts         & 0x00000000000FFFFF;
				uint64_t count_canon     = (counts >> 20) & 0x00000000000FFFFF;
				uint64_t count_non_canon = (counts >> 40) & 0x00000000000FFFFF;
				uint64_t count_both      = count_all - count_canon - count_non_canon;

				// Update general statistics
				only_canonical[count_all][count_canon]++;
				only_non_canonical[count_all][count_non_canon]++;
				both_forms[count_all][count_both]++;

				// Filter or not?
				if(count_all>=minimum_kmer_count) { // pass MAC
					if ((static_cast<double>(count_canon + count_both) >= 
								ceil(minimum_strand_per * static_cast<double>(count_all))) && 
							(static_cast<double>(count_non_canon + count_both) >= 
							 ceil(minimum_strand_per * static_cast<double>(count_all)))) {
						file_out_kmers.write(reinterpret_cast<const char *>(&unique_kmers[i]), sizeof(unique_kmers[i]));
						cnt_pass++;
						shareness[count_all]++; // This statistics only for the used k-mers
					} else {
						file_non_pass_kmers << bits2kmer31(unique_kmers[i], kmer_len) << "\t"
							<< count_all << "\t" << count_canon << "\t" << count_non_canon << "\t"
							<< count_both << endl;
						cnt_no_pass++;
					}
				}	
			}
		}
		cerr << "kmers lower than MAC:\t" << cnt_low_MAC << endl;
		cerr << "passed kmers:\t" << cnt_pass << endl;
		cerr << "passed MAC bot not pass strand filter:\t" << cnt_no_pass << endl;

		file_out_kmers.close();
		file_non_pass_kmers.close();

		/* Output shareness measures */
		ofstream file_shareness(fn_output + ".shareness");
		file_shareness << "kmer appearance\tcount" << endl;
		for(size_t i=0; i<shareness.size(); i++) 
			file_shareness << i << "\t" << shareness[i] << endl;
		file_shareness.close();

		output_uint64_t_matrix(fn_output + ".stats.only_canonical", only_canonical);
		output_uint64_t_matrix(fn_output + ".stats.only_non_canonical", only_non_canonical);
		output_uint64_t_matrix(fn_output + ".stats.both", both_forms);
	} catch (const cxxopts::OptionException& e)
	{
		cerr << "error parsing options: " << e.what() << endl;
		cerr << options.help() << endl;
		exit(1);
	}
	return 0;
}

