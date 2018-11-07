///
///      @file  F_correlate_kmers_to_phenotype.cpp
///     @brief  Loading and correlating phenotype to the presence of k-mers
///
/// We will load a phenotype data and a table of k-mers presence absence.
/// The program will output the interesting k-mers as well as their presence-absence
/// information to output files.
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

#include "kmer_general.h"
#include "kmer_multipleDB.h"

#include <algorithm>    // 
#include <iostream>
#include <iterator>
#include <utility> //std::pair

#include <boost/program_options.hpp>

#include "CTPL/ctpl_stl.h" //thread pool library

using namespace std;
namespace po = boost::program_options;

///
/// @brief  Loading phenotype values from file
/// @param  filename: a path to a file containing phenotypes for each accession
/// @return a phenotype list, pair of two vectors: 1. contain the accessions indices and the second the 
/// phenotype for each accession
///
phenotype_list load_phenotypes_file(const string &filename) {
	phenotype_list p_list;
	
	std::ifstream fin(filename);
	fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // read until the end of the first line
	std::string word;
	double value;
	while(fin >> word) {
		fin >> value;
		p_list.first.push_back(word);
		p_list.second.push_back(value);
	}	
	return p_list;
}

///
/// @brief  Loading phenotypes values from file contatining more than one phenotype
/// @param  filename: a path to a file containing phenotypes for each accession
/// @return vactor of phenotype lists, pair of two vectors: 1. contain the accessions indices and the second the 
/// phenotype for each accession
///

bool is_uint(const std::string& s)
{
	return !s.empty() && std::find_if(s.begin(), 
			s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
}


pair<vector<string>, vector<uint64_phenotype_list>> load_uint64_phenotypes_file(const string &filename) {
	vector<uint64_phenotype_list> p_list;
	vector<string> phenotypes_names;

	std::ifstream 	fin(filename);
	vector<string>  line_tokens;
	string          line, cell;
	size_t line_n(0);
	while(getline(fin, line)) { // Read the file line by line
		stringstream          lineStream(line);
		while(std::getline(lineStream, cell, '\t'))
			line_tokens.push_back(cell);
		// Parsing the file into phenotypes
		if(line_n == 0) {// header
			for(size_t i=1; i<line_tokens.size(); i++)
				phenotypes_names.push_back(line_tokens[i]);
			p_list.resize(phenotypes_names.size()); 
		} else { // data line
			if(line_tokens.size() != (phenotypes_names.size() + 1)) 
				throw std::logic_error("File should have the same number of fields in each row | " +
						filename);
			for(size_t i=0; i<phenotypes_names.size(); i++) {
				p_list[i].first.push_back(line_tokens[0]);
				if(is_uint(line_tokens[i+1]))
					p_list[i].second.push_back(stoul(line_tokens[i+1]));
				else
					throw std::logic_error("Phenotypes should be natural numbers (uints) " + line_tokens[i+1]);
			}
		}
		line_n++;
		line_tokens.resize(0);
	}
	return pair<vector<string>, vector<uint64_phenotype_list>> (phenotypes_names, p_list);
}
// Notice that the names of the accession in each of the phenotypes has to be the same
void write_fam_file(const vector<uint64_phenotype_list> &phenotypes, const string &fn) {
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

void write_fam_file(const uint64_phenotype_list &phenotype, const string &fn) {
	write_fam_file(vector<uint64_phenotype_list>{phenotype}, fn);
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
uint64_phenotype_list intersect_phenotypes_to_present_DBs(const uint64_phenotype_list &pl, const vector<KMC_db_handle> &DB_paths, const bool &must_be_present) {
	uint64_phenotype_list intersect_pl;
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

vector<string> get_DBs_names(const vector<KMC_db_handle> &DBs) {
	vector<string> names;
	for(size_t i=0; i<DBs.size(); i++) {
		names.push_back(DBs[i].name);
	}
	return names;
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
vector<uint64_t> make_list_uint64(const vector<double> &l) {
	vector<uint64_t> res(l.size());
	for(size_t i=0; i<l.size(); i++) {
		if(is_double_uint(l[i]))
			res[i] = (uint64_t)l[i];
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
			("kmers_table",			po::value<string>(),	"Presence/absemce k-mer file")
			("best,n",				po::value<size_t>()->default_value(1000000), 
			 "Number of best k-mers to report")
			("batch_size",			po::value<size_t>()->default_value(10000000), 
			 "Loading only part of the presence absence info to memory")
			("parallel",			po::value<size_t>()->default_value(4), 
			 "Max number of threads to use")
			("kmer_len",			po::value<uint32_t>(), 
			 "Length of the k-mers")
			("maf",			po::value<double>()->default_value(0.05), 
			 "Minor allele frequency")
			("mac",			po::value<size_t>()->default_value(5), 
			 "Minor allele count")
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
		/***************************************************************************************************/
		// Base name for all save files
		double t0 = get_time();
		string fn_base = vm["output_dir"].as<string>() + "/" + vm["base_name"].as<string>();
		size_t heap_size = vm["best"].as<size_t>();	// # best k-mers to report
		size_t batch_size = vm["batch_size"].as<size_t>(); // Load each time ~1/n_steps of k-mers
		ctpl::thread_pool tp(vm["parallel"].as<size_t>());

		uint32_t kmer_length = vm["kmer_len"].as<uint32_t>();
		double maf = vm["maf"].as<double>();
		size_t mac = vm["mac"].as<size_t>();

		// Load DB paths
		vector<KMC_db_handle> DB_paths = read_accession_db_list(vm["paths_file"].as<string>());
		// Loading the phenotype (also include the list of needed accessions)
		pair<vector<string>, vector<uint64_phenotype_list>> phenotypes_info = load_uint64_phenotypes_file(
				vm["phenotype_file"].as<string>());
		size_t phenotypes_n = phenotypes_info.first.size();

		// Intersect phenotypes only to the present DBs (can also check if all must be present)
		for(size_t i=0; i<phenotypes_n; i++)
			phenotypes_info.second[i] = intersect_phenotypes_to_present_DBs(phenotypes_info.second[i], 
					DB_paths, true);
		vector<uint64_phenotype_list> p_list{phenotypes_info.second};

		// Load all accessions data to a combine dataset
		kmer_multipleDB multiDB(
				vm["kmers_table"].as<string>(),
				get_DBs_names(DB_paths),
				p_list[0].first, 
				kmer_length);

		kmer_multipleDB multiDB_step2(
				vm["kmers_table"].as<string>(),
				get_DBs_names(DB_paths),
				p_list[0].first, 
				kmer_length);

		// Create heaps to save all best k-mers & scores
		vector<kmer_heap> k_heap(phenotypes_n, kmer_heap(heap_size)); 
		vector<std::future<void>> tp_results(k_heap.size());
		/****************************************************************************************************
		 *	Load all k-mers and presence/absence information and correlate with phenotypes
		 ***************************************************************************************************/
		size_t min_count = ceil(double(p_list[0].first.size())*maf); // MAF of 5% - maybe should make this parameter external
		if(min_count < mac)
			min_count = mac;

		cerr << "Min count to associate = " << min_count << endl;
		size_t batch_index = 0;
		while(multiDB.load_kmers(batch_size)) { 
			cerr << "Associating k-mers, part: " <<  batch_index << endl; 
			for(size_t j=0; j<(phenotypes_n); j++) { // Check association for each sample 
				tp_results[j] = tp.push([&multiDB,&k_heap,&p_list,j,min_count](int){
						multiDB.add_kmers_to_heap(k_heap[j], p_list[j].second, min_count);});
			}
			for(size_t j=0; j<(phenotypes_n); j++) {
				tp_results[j].get();
				cerr << ".";
				cerr.flush();
			}
			cerr << endl;
			batch_index++;
		}

		cerr << "[TIMING]\tAFTER LOADING\t" << ((get_time() - t0)/60) << endl;
		// close DBs, and release space ??
		/****************************************************************************************************
		 *	save best k-mers to files and create set of k-mers from heaps (and clear heaps)	
		 ***************************************************************************************************/
		vector<kmers_output_list> best_kmers;	
		cerr << "[TIMING]\tBEFORE SAVING\t" << ((get_time() - t0)/60) << endl;
		for(size_t j=0; j<(phenotypes_n); j++) {
			string fn_kmers = fn_base + "." + std::to_string(j) + ".best_kmers";
			k_heap[j].output_to_file_with_scores(fn_kmers + ".scores");

			// Create set of best k-mers
			best_kmers.push_back(k_heap[j].get_kmers_for_output(kmer_length));
			k_heap[j].empty_heap(); // clear heap
		}
		cerr << "[TIMING]\tAFTER SAVING\t" << ((get_time() - t0)/60) << endl;

		/****************************************************************************************************
		 *	Reload all k-mers and out to plink files only the best k-mers found in the previous stage
		 ***************************************************************************************************/
		// Open all bed & bim plink files	
		vector<bedbim_handle> plink_output;
		for(size_t j=0; j<(phenotypes_n); j++) {
			plink_output.emplace_back(fn_base + "." + std::to_string(j) + "." + phenotypes_info.first[j] ) ;
			write_fam_file(p_list[j], fn_base + "." + std::to_string(j) + "." + phenotypes_info.first[j] 
					+ ".fam");
		}
		cerr << "[TIMING]\tOPENED PLINK FILES\t" << ((get_time() - t0)/60) << endl;

		// Reload k-mers to create plink bed/bim files
		batch_index = 0;
		while(multiDB_step2.load_kmers(batch_size)) { 
			cerr << "Saving k-mers, part: " <<  batch_index << endl; 
			for(size_t j=0; j<(phenotypes_n); j++) { // Check association for each samples  
				best_kmers[j].next_index = 
					multiDB_step2.output_plink_bed_file(plink_output[j], best_kmers[j].list, best_kmers[j].next_index); // parallel ?
				cerr << ".";
				cerr.flush();	
			}
			cerr << endl;
			batch_index++;
		}
		cerr << "[TIMING]\tWROTE PLINK FILES\t" << ((get_time() - t0)/60) << endl;
		// close DBs, release space ??
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


