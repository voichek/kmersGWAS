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
#include "kmers_multiple_databases.h"
#include "best_associations_heap.h"

#include <sys/sysinfo.h> // To monitor memory usage
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
pair<vector<string>, vector<PhenotypeList>> load_phenotypes_file(const string &filename) {
	vector<PhenotypeList> p_list;
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
				p_list[i].second.push_back(stof(line_tokens[i+1]));
			}
		}
		line_n++;
		line_tokens.resize(0);
	}
	return pair<vector<string>, vector<PhenotypeList>> (phenotypes_names, p_list);
}
// Notice that the names of the accession in each of the phenotypes has to be the same
void write_fam_file(const vector<PhenotypeList> &phenotypes, const string &fn) {
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

void write_fam_file(const PhenotypeList &phenotype, const string &fn) {
	write_fam_file(vector<PhenotypeList>{phenotype}, fn);
}

size_t get_index_DB(const string &name, const vector<KMCDataBaseHandle> &DBs) {
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

PhenotypeList intersect_phenotypes_to_present_DBs(const PhenotypeList &pl, const vector<KMCDataBaseHandle> &DB_paths, const bool &must_be_present) {
	PhenotypeList intersect_pl;
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

vector<string> get_DBs_paths(const vector<string> &names, const vector<KMCDataBaseHandle> &DBs) {
	vector<string> paths;
	for(size_t i=0; i<names.size(); i++) {
		size_t index = get_index_DB(names[i], DBs);
		if(index == (~0u))
			throw std::logic_error("Couldn't find DB " + names[i]);
		paths.push_back(DBs[index].dir_path);
	}
	return paths;
}

vector<string> get_DBs_names(const vector<KMCDataBaseHandle> &DBs) {
	vector<string> names;
	for(size_t i=0; i<DBs.size(); i++) {
		names.push_back(DBs[i].name);
	}
	return names;
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
			("k_mers_scores", "output the best k_mers scores in binary format")
			("debug_option_batches_to_run",			po::value<uint64_t>()->default_value(NULL_KEY), 
			 "Change to run only on part of the table")
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
		string fn_base = vm["output_dir"].as<string>() + "/" + vm["base_name"].as<string>();
		size_t heap_size = vm["best"].as<size_t>();	// # best k-mers to report
		size_t batch_size = vm["batch_size"].as<size_t>(); // Load each time ~1/n_steps of k-mers
		ctpl::thread_pool tp(vm["parallel"].as<size_t>());

		uint32_t kmer_length = vm["kmer_len"].as<uint32_t>();
		double maf = vm["maf"].as<double>();
		size_t mac = vm["mac"].as<size_t>();

		uint64_t debug_option_batches_to_run = vm["debug_option_batches_to_run"].as<uint64_t>();
		// Load DB paths
		vector<KMCDataBaseHandle> DB_paths = read_accession_db_list(vm["paths_file"].as<string>());
		// Loading the phenotype (also include the list of needed accessions)
		pair<vector<string>, vector<PhenotypeList>> phenotypes_info = load_phenotypes_file(
				vm["phenotype_file"].as<string>());
		size_t phenotypes_n = phenotypes_info.first.size();

		// Intersect phenotypes only to the present DBs (can also check if all must be present)
		for(size_t i=0; i<phenotypes_n; i++)
			phenotypes_info.second[i] = intersect_phenotypes_to_present_DBs(phenotypes_info.second[i], 
					DB_paths, true);
		vector<PhenotypeList> p_list{phenotypes_info.second};

		// Load all accessions data to a combine dataset
		MultipleKmersDataBases multiDB(
				vm["kmers_table"].as<string>(),
				get_DBs_names(DB_paths),
				p_list[0].first, 
				kmer_length);

		MultipleKmersDataBases multiDB_step2(
				vm["kmers_table"].as<string>(),
				get_DBs_names(DB_paths),
				p_list[0].first, 
				kmer_length);

		// Create heaps to save all best k-mers & scores
		vector<BestAssociationsHeap> k_heap(phenotypes_n, BestAssociationsHeap(heap_size)); 
		vector<std::future<void>> tp_results(k_heap.size());
		/****************************************************************************************************
		 *	Load all k-mers and presence/absence information and correlate with phenotypes
		 ***************************************************************************************************/
		size_t min_count = ceil(double(p_list[0].first.size())*maf); // MAF of 5% - maybe should make this parameter external
		if(min_count < mac)
			min_count = mac;

		double t0,t1;
		/* Compute time taken */
		cerr << "Min count to associate = " << min_count << endl;
		cerr << "Used RAM:\t" << get_mem_used_by_process() << endl;
		size_t batch_index = 0;
		t0 = get_time();
		while(multiDB.load_kmers(batch_size) && (batch_index < debug_option_batches_to_run)) { 
			t1 = get_time();
			cerr << "Associating k-mers, part: " <<  batch_index << "\tt(min)=" << (double)(t1-t0)/(60.) << endl; 
			cerr << "Used RAM:\t" << get_mem_used_by_process() << endl;
			t0 = get_time();
			for(size_t j=0; j<(phenotypes_n); j++) { // Check association for each sample 
				tp_results[j] = tp.push([&multiDB,&k_heap,&p_list,j,min_count](int){
						multiDB.add_kmers_to_heap(k_heap[j], p_list[j].second, min_count);});
			}
			for(size_t j=0; j<(phenotypes_n); j++) {
				tp_results[j].get();
				cerr << ".";
				cerr.flush();
			}
			t1 = get_time();
			cerr << "\tt(min)="<< (double)(t1-t0)/(60.) <<endl;
			t0 = get_time();
			batch_index++;
		}

		cerr << "Used RAM:\t" << get_mem_used_by_process() << endl;
		// close DBs, and release space ??
		/****************************************************************************************************
		 *	save best k-mers to files and create set of k-mers from heaps (and clear heaps)	
		 ***************************************************************************************************/
		vector<kmers_output_list> best_kmers;	
		for(size_t j=0; j<(phenotypes_n); j++) { // For some reason this can be improved by parallelization
			if (vm.count("k_mers_scores")) {
				string fn_kmers = fn_base + "." + std::to_string(j) + ".best_kmers";
				cerr << "output: " << fn_kmers << "\t\t";
				cerr << "Used RAM:\t" << get_mem_used_by_process() ;
				k_heap[j].output_to_file_with_scores(fn_kmers + ".scores");
				cerr << "\t" << get_mem_used_by_process();
			}
			// Create set of best k-mers
			best_kmers.push_back(k_heap[j].get_kmers_for_output(kmer_length));
			cerr << "\t" << get_mem_used_by_process();
			k_heap[j].empty_heap(); // clear heap
			cerr << "\t" << get_mem_used_by_process() << endl;
		}
		cerr << "Used RAM:\t" << get_mem_used_by_process() << endl;

		/****************************************************************************************************
		 *	Reload all k-mers and out to plink files only the best k-mers found in the previous stage
		 ***************************************************************************************************/
		// Open all bed & bim plink files	
		vector<BedBimFilesHandle> plink_output;
		for(size_t j=0; j<(phenotypes_n); j++) {
			plink_output.emplace_back(fn_base + "." + std::to_string(j) + "." + phenotypes_info.first[j] ) ;
			write_fam_file(p_list[j], fn_base + "." + std::to_string(j) + "." + phenotypes_info.first[j] 
					+ ".fam");
		}
		cerr << "Re-loading k-mers" << endl;
		cerr << "Used RAM:\t" << get_mem_used_by_process() << endl;
		// Reload k-mers to create plink bed/bim files
		batch_index = 0;
		while(multiDB_step2.load_kmers(batch_size) && (batch_index<debug_option_batches_to_run)) {
			cerr << "Saving k-mers, part: " <<  batch_index << endl; 
			cerr << "Used RAM:\t" << get_mem_used_by_process() << endl;
			for(size_t j=0; j<(phenotypes_n); j++) { // Check association for each samples  
				best_kmers[j].next_index = 
					multiDB_step2.output_plink_bed_file(plink_output[j], 
							best_kmers[j].list, best_kmers[j].next_index); // parallel ?
				cerr << ".";
				cerr.flush();	
			}
			cerr << endl;
			batch_index++;
		}
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


