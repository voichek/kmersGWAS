///
///      @file  kmers_add_strand_information.cpp
///     @brief  Take the KMC k-mers counts with and without cononization and add to the k-mers
///     on which strands there were found
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  03/10/19
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2019, Yoav Voichek
///
///This source code is released for free distribution under the terms of the
///GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmer_general.h"

#include <bitset>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <cxxopts/include/cxxopts.hpp>

using namespace std;

inline tuple<uint64_t, uint64_t> is_canonized_kmer_representation_flag(const uint64_t& k, const uint32& k_len) {
	uint64_t k_rc = kmer_reverse_complement(k, k_len); 
	if(k<k_rc) 
		return make_tuple(k,    0x4000000000000000);
	else
		return make_tuple(k_rc, 0x8000000000000000);
}

int main(int argc, char* argv[]) {
	cxxopts::Options options("kmers_add_strand_information", "Combines the two KMC DB (canonized and not-canonized) and output the list of k-mers in binary format");
	try
	{
		options.add_options()
			("c,canonized_KMC", "prefix of KMC DB files, run with canonization of k-mers", cxxopts::value<string>())
			("n,non_canonized_KMC", "prefix of KMC DB files, run without canonization of k-mers", cxxopts::value<string>())
			("k,kmers_len", "length of k-mers", cxxopts::value<uint32>())
			("o,output", "output file name", cxxopts::value<string>())
			("help", "print help")
			;
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cerr << options.help() << endl;
			exit(0);
		}
		vector<string> required_parametrs({"canonized_KMC", "non_canonized_KMC", "kmers_len", "output"});
		for(size_t i=0; i<required_parametrs.size(); i++) {
			if(result.count(required_parametrs[i]) == 0) {
				cerr << required_parametrs[i] << " is a required parameter" << endl;
				cerr << options.help() << endl;
				exit(1);
			}
		}
		// Load parameters
		string fn_prefix_KMC_canon(result["canonized_KMC"].as<string>());
		string fn_prefix_KMC_non_canon(result["non_canonized_KMC"].as<string>());
		string fn_output(result["output"].as<string>());
		uint32 kmer_length(result["kmers_len"].as<uint32>());

		// Check if all input files exist
		vector<string> required_files({fn_prefix_KMC_canon + ".kmc_pre", fn_prefix_KMC_canon + ".kmc_suf",
				fn_prefix_KMC_non_canon + ".kmc_pre", fn_prefix_KMC_non_canon + ".kmc_suf"});
		for(size_t i=0; i<required_files.size(); i++) {
			if(!is_file_exist(required_files[i])) {
				cerr << "Couldn't find file: " << required_files[i] << endl;
				exit(1);
			}
		}
		if((kmer_length > 31) || (kmer_length <10)) {
			cerr << "kmer length has to be between 10-31" << endl;
			exit(1);
		}
		/****************************************************************************************************/
		/* END of parsing and checking input parameters */
		/****************************************************************************************************/
		// Go over the regular KMC database kmers and add to hashtable as keys with value 0
		KmerUint64Hash kmers(1000000); // Can get the initial hash_table_size from user
		kmers.set_empty_key(NULL_KEY);


		CKmerUpTo31bpAPI kmer_obj(kmer_length); uint kmer_counter;
		CKMCFile kmer_db_canon;
		kmer_db_canon.OpenForListing(fn_prefix_KMC_canon);
		while(kmer_db_canon.ReadNextKmer(kmer_obj, kmer_counter)) 
			kmers.insert(KmerUint64Hash::value_type(kmer_obj.to_uint(),0));
		cout << "Canonized kmers:\t" << kmers.size() << endl;

		// Go over the non-canon' KMC database save flags (0x1 / 0x2) if apeared as canon/non-canon
		CKMCFile kmer_db_non_canon;
		kmer_db_non_canon.OpenForListing(fn_prefix_KMC_non_canon);
		KmerUint64Hash::iterator it_hash;
		size_t counter_all(0), counter_found(0);
		while(kmer_db_non_canon.ReadNextKmer(kmer_obj, kmer_counter)) {
			uint64_t k, flag;
			tie(k, flag)  = is_canonized_kmer_representation_flag(kmer_obj.to_uint(), kmer_length);
			it_hash = kmers.find(k);
			if(it_hash != kmers.end()) {
				it_hash->second |= flag;
				counter_found++;
			}
			counter_all++;
		}

		cout << "Non-canon kmers:\t" << counter_all << endl;
		cout << "Non-canon kmers found:\t" << counter_found << endl;
		// output to file all kmers with the flags in the most significant bits (as k-mers is of len < 31
		// there is at least two free bit
		vector<size_t> flag_counter(4,0);
		vector<uint64_t> kmers_to_save;
		for(auto it : kmers) {
			flag_counter.at(it.second >> 62)++;
			kmers_to_save.push_back(it.first | it.second); // add flags in 2-MSBs
		}
		for(size_t i=0; i<4; i++)
			cout << "flag\t" << i << "\tcount is\t" << flag_counter[i] << endl;
		

		if(flag_counter[0] != 0) {
			// When counting k-mers without canonization all k-mers should be counted without a threshold on minimal appearences
			cerr << "Error: flag 00 should be equal to zero." << endl;
			cerr << "This is likely due to running the KMC non-canonized with -ci not 1" << endl;
			exit(1);
		}
			

		sort(kmers_to_save.begin(), kmers_to_save.end(),
				[](const uint64_t & a, const uint64_t & b) -> bool
				{return (a&0x3FFFFFFFFFFFFFFF) < (b&0x3FFFFFFFFFFFFFFF);});
		cout << "kmers to save:\t" << kmers_to_save.size() << endl;

		ofstream fout(fn_output, ios::binary);
		for(size_t i=0; i<kmers_to_save.size(); i++)
			fout.write(reinterpret_cast<const char *>(&kmers_to_save[i]), sizeof(kmers_to_save[i]));
		fout.close();
	} catch (const cxxopts::OptionException& e)
	{
		cerr << "error parsing options: " << e.what() << endl;
		cerr << options.help() << endl;
		exit(1);
	}
	return 0;
}

