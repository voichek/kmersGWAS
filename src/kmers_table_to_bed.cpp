///
///      @file  kmers_table_to_bed.cpp
///     @brief  Take a kmers table and convert it to bed/bim file(s)
///
/// 	In some case we would like to run the full kmers presence/absence information in the exact
/// 	LMM method of choice (e.g. GEMMA). For this purpose we need to convert the kmers table to 
/// 	bed/bim/fam format, which is the purpose of this code. To save space, we will filter the 
/// 	kmers to have the wanted minor allele frequency / count (MAF, MAC) before outputing them.
/// 	Moreover, we will also filter kmers which have the same exact presence/absence pattern.
/// 	This will be done to reduce the computational load on the next step (e.g. GEMMA).
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  12/23/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include <iostream>
#include <string>
#include <vector>

#include <cxxopts/include/cxxopts.hpp>

#include "kmer_general.h"
#include "kmers_multiple_databases.h"


using namespace std;

int main(int argc, char* argv[]) {
	cxxopts::Options options("kmers_table_to_bed", "Convert k-mers table to PLINK binary format");
	try
	{
		options.add_options()
			("t,kmers_table", "k-mers table path", cxxopts::value<string>())
			("k,kmers_len", "length of k-mers", cxxopts::value<size_t>())
			("p,phentype_file", "phenotype file, condense output only to individuals with a phenotype", cxxopts::value<string>())
			("maf", "minor allele frequency", cxxopts::value<double>())
			("mac", "minor allele count", cxxopts::value<size_t>())
			("b,batch_size", "maximal number of variants in each PLINK bed file (seperate to many file  if needed)", cxxopts::value<size_t>())
			("o,output", "prefix for output files", cxxopts::value<string>())
			("u,unique_patterns", "output only unique presence/absence patterns", cxxopts::value<bool>()->default_value("false"))
			("help", "print help")
			;
		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cerr << options.help() << endl;
			exit(0);
		}
		vector<string> required_parametrs({"kmers_table", "kmers_len", "phentype_file", "maf", "mac", "batch_size", "output"});
		for(size_t i=0; i<required_parametrs.size(); i++) {
			if(result.count(required_parametrs[i]) == 0) {
				cerr << required_parametrs[i] << " is a required parameter" << endl;
				cerr << options.help() << endl;
				exit(1);
			}
		}
		// Load parameters
		string fn_kmers_table(result["kmers_table"].as<string>());
		string fn_phenotypes(result["phentype_file"].as<string>());
		double MAF(result["maf"].as<double>());
		size_t MAC(result["mac"].as<size_t>());
		size_t max_batch_size(result["batch_size"].as<size_t>());
		size_t kmer_len(result["kmers_len"].as<size_t>());
		string output_base(result["output"].as<string>());
		
		bool unique_patterns(result["unique_patterns"].as<bool>());
		// Check if all input files exist
		vector<string> required_files({fn_kmers_table+".names",fn_kmers_table+".table",fn_phenotypes});
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
		
		// load phenotypes (needed for current DBs)
		pair<vector<string>, vector<PhenotypeList>> phenotypes_file_info = load_phenotypes_file(fn_phenotypes);
		cerr << "using " << phenotypes_file_info.first[0] << endl;
		PhenotypeList pheno_info = intersect_phenotypes_to_present_DBs(phenotypes_file_info.second[0], fn_kmers_table, true);

		// calculate effective MAC (if MAF indicate lower number)
		size_t min_count = ceil(double(pheno_info.first.size())*MAF); // MAF of 5% - maybe should make this parameter external
		if(min_count < MAC)
			min_count = MAC;

		// create kmers-multiple dbs file
		MultipleKmersDataBases multiDB(fn_kmers_table, pheno_info.first, kmer_len);

		// create presence/absence pattern counter (initialize with 100M)
		KmersSet pa_patterns_counter(1000*1000*100); 
		pa_patterns_counter.set_empty_key(NULL_KEY); // need to define empty value for google dense hash table

		size_t batch_index(0);
		cerr << "loading.... " << endl;
		while(multiDB.load_kmers( max_batch_size, min_count)) { // Running over the kmers table in batches
			cerr << "Batch:\t" << batch_index+1 << endl;
			// output current batch
			string fn_bedbim = output_base + "." + std::to_string(batch_index);
			if(unique_patterns) {
				multiDB.output_plink_bed_file_unique_presence_absence_patterns(fn_bedbim, pa_patterns_counter);
			} else {
				multiDB.output_plink_bed_file(fn_bedbim);
			}
			// output fam file (the same in all repeats)
			ofstream fout(fn_bedbim + ".fam");
			for(size_t i=0; i<pheno_info.first.size(); i++) 
				fout << pheno_info.first[i] << " " <<
					pheno_info.first[i] << " 0 0 0 " <<
					pheno_info.second[i] << endl;

			fout.close();
			batch_index++; // advance batch counter
		}	
	} catch (const cxxopts::OptionException& e)
	{
		cerr << "error parsing options: " << e.what() << endl;
		cerr << options.help() << endl;
		exit(1);
	}
	return 0;
}
