///
///      @file  F_correlate_kmers_to_phenotype.cpp
///     @brief  Loading and correlating phenotype to the presence of k-mers
///
/// We will load a phenotype data and a list of k-mers DBs and will correlate them.
/// The program will output the interesting k-mers as well as their presence-absence
/// information to output files.
/// The program will have the possibility of looking on only part of the k-mers spectrum.	
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  07/17/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///


/**
 * @brief 19/07/2018 - How to build for prototype
 * 
 * Initialy I am going to run all in one program no threads nor parallelization
 * Also all phenotypes and accession will be the same and ordered
 * all the needed complications will come in the next goes
 */

#include <utility> //std::pair
#include "kmer_general.h"
#include "kmer_DB.h"
#include "kmer_multipleDB.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
#include <algorithm>    // std::random_shuffle
#include <cstdlib>      // std::rand, std::srand

#include "CTPL/ctpl_stl.h" //thread pool library
using namespace std;



///
/// @brief  Loading phenotype values from file
/// @param  filename: a path to a file containing phenotypes for each accession
/// @return a phenotype list, pair of two vectors: 1. contain the accessions indices and the second the 
/// phenotype for each accession
///
phenotype_list load_phenotypes_file(const string &filename) {
	phenotype_list p_list;
	
	std::ifstream fin(filename);
	fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::string word;
	double value;
	while(fin >> word) {
		fin >> value;
		p_list.first.push_back(word);
		p_list.second.push_back(value);
	}	
	return p_list;
}

phenotype_list randomize_phenotype(const phenotype_list &orig) {
	phenotype_list permute_values(orig); // copy the original values;
	std::random_shuffle(permute_values.second.begin(), permute_values.second.end());
	return permute_values;
}

// Notice that the names of the accession in each of the phenotypes has to be the same
void write_fam_file(const vector<phenotype_list> &phenotypes, const string &fn) {
	ofstream f(fn, ios::out);
	for(size_t i=0; i<phenotypes[0].first.size(); i++) {
		f << phenotypes[0].first[i] << " " << phenotypes[0].first[i] << " 0 0 0";
		for(size_t j=0; j<phenotypes.size(); j++) {
			if(phenotypes[j].first[i] != phenotypes[0].first[i]) { // check phenotypes have the same order
				throw std::logic_error("phenotypes should have the same order " + 
						phenotypes[j].first[i] + "!=" + phenotypes[0].first[i]);
			}
			f << " " << phenotypes[j].second[i];
		}
		f << endl;
	}
	f.close();
}

void write_fam_file(const phenotype_list &phenotype, const string &fn) {
	write_fam_file(vector<phenotype_list>{phenotype}, fn);
}

size_t get_index_DB(const string &name, const vector<KMC_db_handle> &DBs) {
	size_t index = (~0u);
	for(size_t j=0; j<DBs.size(); j++) {
		if(DBs[j].name == name) {
			if(index != (~0u))
				throw std::logic_error("Two DBs with the same name! " + name);
			index = j;
		}
	}
	return index;
}
phenotype_list intersect_phenotypes_to_present_DBs(const phenotype_list &pl, const vector<KMC_db_handle> &DB_paths, const bool &must_be_present) {
	phenotype_list intersect_pl;
	for(size_t i=0; i<pl.first.size(); i++) {
		size_t index = get_index_DB(pl.first[i], DB_paths);
		if(index == (~0u)) { // nothing found
			if(must_be_present)
				throw std::logic_error("Couldn't find path for DB: " + pl.first[i]);
		} else {
			intersect_pl.first.push_back(pl.first[i]);
			intersect_pl.second.push_back(pl.second[i]);
		}
	}
	return intersect_pl;	
}

vector<string> get_DBs_paths(const vector<string> &names, const vector<KMC_db_handle> &DBs) {
	vector<string> paths;
	for(size_t i=0; i<names.size(); i++) {
		size_t index = get_index_DB(names[i], DBs);
		if(index == (~0u))
			throw std::logic_error("Couldn't find DB " + names[i]);
		paths.push_back(DBs[index].dir_path);
	}
	return paths;
}

bool is_double_uint(const double &x) {
	if(x<0)
		return false;
	if(trunc(x) == x)
		return true;
	else
		return false;
}

// This is
vector<uint64> make_list_uint64(const vector<double> &l) {
	vector<uint64> res(l.size());
	for(size_t i=0; i<l.size(); i++) {
		if(is_double_uint(l[i]))
			res[i] = (uint64)l[i];
		else
			throw std::logic_error("Got a value " + to_string(l[i]) + " need to get unsigned integers");
	}

	return res;

}

int main(int argc, char* argv[])
{
	/*******************************************************************************************************/
	/* Loading the user defined parameters */
	try {
		/* Define the input params */
		po::options_description desc("Allowed options");
		desc.add_options()
			("help", "produce help message")
			("phenotype_file,p",	po::value<string>(),	"phenotype file name")
			("base_name,b",			po::value<string>(),	"base name to use for all files")
			("output_dir,o",		po::value<string>()->default_value("."), 
			 "where to save output files")
			("paths_file",			po::value<string>(),	"file contatining a list of path for every DB")
			("kmers_file",			po::value<string>(),	"name of k-mer file inside every sub-directory")
			("permutations",		po::value<size_t>()->default_value(0), 
			 "How many permutation to do for phenotypes")
			("best,n",				po::value<size_t>()->default_value(1000000), 
			 "Number of best k-mers to report")
			("kmers_parts",			po::value<size_t>()->default_value(500), 
			 "Loading only 1/parts of the k-mers to the memory each time")
			("parallel",			po::value<size_t>()->default_value(4), 
			 "Max number of threads to use")
			;

		/* parse the command line */
		po::variables_map vm;        
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);    

		/* Output help funciton */
		if (vm.count("help")) {
			cout << desc << "\n";
			return 0;
		}

		if (vm.count("phenotype_file")) {
			cerr << "phenotype file: " 
				<< vm["phenotype_file"].as<string>() << "\n";
		} else {
			cout << "Need to specify the phenotype file\n";
			return 1;
		}

		if (vm.count("base_name")) {
			cerr << "base name: " 
				<< vm["base_name"].as<string>() << "\n";
		} else {
			cout << "Need to specify the base name\n";
			return 1;
		}
		//		std::srand ( unsigned ( std::time(0) ) ); 
		std::srand(123456789); // to have determenistic results - important in the premutations
		/***************************************************************************************************/
		cerr << "test" << endl;	

		double t0 = get_time();
		// Base name for all save files
		string fn_base = vm["output_dir"].as<string>() + "/" + vm["base_name"].as<string>();

		size_t perm_phenotype = vm["permutations"].as<size_t>(); // #permutation on phenotypes
		size_t heap_size = vm["best"].as<size_t>();	// # best k-mers to report
		size_t n_steps = vm["kmers_parts"].as<size_t>(); // Load each time ~1/n_steps of k-mers
		ctpl::thread_pool tp(vm["parallel"].as<size_t>());

		// Load DB paths
		vector<KMC_db_handle> DB_paths = read_accession_db_list(vm["paths_file"].as<string>());
		// Loading the phenotype (also include the list of needed accessions)
		vector<phenotype_list> p_list{load_phenotypes_file(vm["phenotype_file"].as<string>())};
		// Intersect phenotypes only to the present DBs (can also check if all must be present)
		p_list[0] = intersect_phenotypes_to_present_DBs(p_list[0], DB_paths, false);

		for(size_t i=0; i<perm_phenotype; i++) // Add permuted phenotypes
			p_list.push_back(randomize_phenotype(p_list[0]));
		write_fam_file(p_list, fn_base + ".fam"); // save all the permutation to a fam file

		// Load all accessions data to a combine dataset
		//		kmer_multipleDB multiDB(vm["DBs_path"].as<string>(), p_list[0].first, vm["kmers_file"].as<string>());    
		kmer_multipleDB multiDB(get_DBs_paths(p_list[0].first, DB_paths), p_list[0].first, vm["kmers_file"].as<string>());
		kmer_multipleDB multiDB_step2(get_DBs_paths(p_list[0].first, DB_paths), p_list[0].first, vm["kmers_file"].as<string>());

		// Create heaps to save all best k-mers & scores
		vector<kmer_heap> k_heap(perm_phenotype+1, kmer_heap(heap_size)); 
		vector<std::future<void>> tp_results(k_heap.size());
		/****************************************************************************************************
		 *	Load all k-mers and presence/absence information and correlate with phenotypes
		 ***************************************************************************************************/
		cerr << "[TIMING]\tBEFORE LOADING\t" << ((get_time() - t0)/60) << endl;
		size_t min_count = (p_list[0].first.size()+20-1)/20;
		cerr << "min_count = " << min_count << endl;
		for(uint64 i=1; i<=n_steps; i++) { 
			multiDB.load_kmers(i, n_steps); 
			for(size_t j=0; j<(perm_phenotype + 1); j++) { // Check association for each permutation
//				multiDB.add_kmers_to_heap(k_heap[j],  // parallel?
//						p_list[j].second, p_list[j].first); // first is same for all j
//				multiDB.add_kmers_to_heap(k_heap[j],  // parallel?
//						make_list_uint64(p_list[j].second)); // first is same for all j
				tp_results[j] = tp.push([&multiDB,&k_heap,&p_list,j,min_count](int){
						multiDB.add_kmers_to_heap(k_heap[j], make_list_uint64(p_list[j].second),min_count);});
			}
			for(size_t j=0; j<(perm_phenotype+1); j++) {
				tp_results[j].get();
				cerr << i << "\t" <<  j << "\t"; k_heap[j].plot_stat(); // heap status
			}
		}
		cerr << "[TIMING]\tAFTER LOADING\t" << ((get_time() - t0)/60) << endl;
		// close DBs, and release space
		multiDB.clear_hashtable();
		/****************************************************************************************************
		 *	save best k-mers to files and create set of k-mers from heaps (and clear heaps)	
		 ***************************************************************************************************/
		vector<kmer_set> best_kmers;	
		cerr << "[TIMING]\tBEFORE SAVING\t" << ((get_time() - t0)/60) << endl;
		for(size_t j=0; j<(perm_phenotype + 1); j++) {
			string fn_kmers = fn_base + "." + std::to_string(j) + ".best_kmers";
			k_heap[j].output_to_file_with_scores(fn_kmers + ".scores");

			// Create set of best k-mers
			best_kmers.push_back(k_heap[j].get_kmer_set());
			k_heap[j].empty_heap(); // clear heap
		}
		cerr << "[TIMING]\tAFTER SAVING\t" << ((get_time() - t0)/60) << endl;

		/****************************************************************************************************
		 *	Reload all k-mers and out to plink files only the best k-mers found in the previous stage
		 ***************************************************************************************************/
		// Open all bed & bim plink files	
		vector<bedbim_handle> plink_output;
		for(size_t j=0; j<(perm_phenotype+1); j++) {
			plink_output.emplace_back(fn_base + "." + std::to_string(j)) ;
			write_fam_file(p_list[j], fn_base + "." + std::to_string(j) + ".fam");
		}
		cerr << "[TIMING]\tOPENED PLINK FILES\t" << ((get_time() - t0)/60) << endl;

		// Reload k-mers to create plink bed/bim files
		for(uint64 i=1; i<=n_steps; i++) { // Notice here index from 1 and not 0 (might worth changing) 
			multiDB_step2.load_kmers(i, n_steps); 
			for(size_t j=0; j<(perm_phenotype + 1); j++) { // Check association for each permutation
				multiDB_step2.output_plink_bed_file(plink_output[j], best_kmers[j]); // parallel ?
				cerr << i << "\t" <<  j << endl; 
			}
		}
		cerr << "[TIMING]\tWROTE PLINK FILES\t" << ((get_time() - t0)/60) << endl;
		multiDB_step2.clear_hashtable(); // close DBs, release space
		// close all bed & bim plink files & (specific fam?)
		for(size_t j=0; j<plink_output.size(); j++)
			plink_output[j].close();

	}
	catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	}
	catch(...) {
		cerr << "Exception of unknown type!\n";
	}
	return 0;
}


