///
///      @file  YV_correlate_kmers_to_phenotype.cpp
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
using namespace std;

typedef pair <vector<string>, vector<double> > phenotype_list;


///
/// @brief  Loading phenotype values from file
/// @param  filename: a path to a file containing phenotypes for each accession
/// @return a phenotype list, pair of two vectors: 1. contain the accessions indices and the second the 
/// phenotype for each accession
///
phenotype_list load_phenotypes_file(string filename) {
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

int main(int argc, char* argv[])
{
	/*******************************************************************************************************/
	/* Loading the user defined parameters */
	try {
		/* Define the input params */
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("phenotype_file,p", po::value<string>(), "path to the phenotype file")
            ("DBs_path,d", po::value<string>()->default_value("/tmp/global2/yvoichek/kmer_counts/"), 
			 "path of the k-mer DBS")
            ("kmers_file", po::value<string>()->default_value("order_kmers_appear_more_than_once"), 
			 "path of the k-mer DBS")
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
                 << vm["phenotype_file"].as<string>() << ".\n";
        } else {
            cout << "Need to specify the phenotype file\n";
			return 1;
        }
		
		/***************************************************************************************************/
		// 1. Loading the phenotype (also include the list of needed accessions)
		phenotype_list p_list = load_phenotypes_file(vm["phenotype_file"].as<string>());

		// 2. Load all accessions data to a combine dataset
		kmer_multipleDB multiDB(vm["DBs_path"].as<string>(), p_list.first, vm["kmers_file"].as<string>());    
		
		// 3. test
		double t0, t1;
		for(uint64 i=1; i<=100; i++) {
			t0 = get_time();
			multiDB.load_kmers(i, 100);
			t1 = get_time();
			cout << i << "\t" << (t1-t0)/60 << endl;
		}
//		multiDB.plot_textual_hash_map();
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
