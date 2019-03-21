
#include "kmer_general.h"
#include "kmers_multiple_databases.h"

#include <sys/sysinfo.h> // To monitor memory usage
#include <algorithm>    // 
#include <iostream>
#include <iterator>
#include <utility> //std::pair

#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;


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
			("paths_file",			po::value<string>(),	"file contatining a list of path for every DB")
			("kmers_table",			po::value<string>(),	"Presence/absemce k-mer file")
			("batch_size",			po::value<size_t>()->default_value(10000000), 
			 "Loading only part of the presence absence info to memory")
			("kmer_len",			po::value<uint32_t>(), 
			 "Length of the k-mers")
			("maf",			po::value<double>()->default_value(0.05), 
			 "Minor allele frequency")
			("mac",			po::value<size_t>()->default_value(5), 
			 "Minor allele count")
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
		/***************************************************************************************************/
		// Base name for all save files
		size_t batch_size = vm["batch_size"].as<size_t>(); // Load each time ~1/n_steps of k-mers

		uint32_t kmer_length = vm["kmer_len"].as<uint32_t>();
		double maf = vm["maf"].as<double>();
		size_t mac = vm["mac"].as<size_t>();

		uint64_t debug_option_batches_to_run = vm["debug_option_batches_to_run"].as<uint64_t>();


		// Loading the phenotype (also include the list of needed accessions)
		pair<vector<string>, vector<PhenotypeList>> phenotypes_info = load_phenotypes_file(
				vm["phenotype_file"].as<string>());

		// Intersect phenotypes only to the present DBs (can also check if all must be present)
		phenotypes_info.second[0] = intersect_phenotypes_to_present_DBs(phenotypes_info.second[0],
			vm["kmers_table"].as<string>(), true);
		vector<PhenotypeList> p_list{phenotypes_info.second};

		// Load all accessions data to a combine dataset
		MultipleKmersDataBases multiDB(
				vm["kmers_table"].as<string>(),
				p_list[0].first, 
				kmer_length);

		size_t n_acc = p_list[0].first.size();
		size_t min_count = ceil(double(p_list[0].first.size())*maf); // MAF of 5% - maybe should make this parameter external
		if(min_count < mac)
			min_count = mac;

		double t0,t1;
		/* Compute time taken */
		cerr << "Min count to associate = " << min_count << endl;
		size_t batch_index = 0;

		size_t M(0);
		vector<vector<double> > R(n_acc, vector<double>(n_acc, 0));

		t0 = get_time();
		while(multiDB.load_kmers(batch_size, min_count) && (batch_index < debug_option_batches_to_run)) { 
			t1 = get_time();
			cerr << "calc pre-gamma params part: " <<  batch_index << "\tt(min)=" << (double)(t1-t0)/(60.) << endl; 
			t0 = get_time();

			multiDB.update_gamma_precalculations(R, M);
			t1 = get_time();
			cerr << "\tt(min)="<< (double)(t1-t0)/(60.) << "\tM="<<M <<endl;
			t0 = get_time();
			batch_index++;
		}
		// copy lower triangle to upper one in R
		for(size_t i=0; i<n_acc; i++) {
			for(size_t j=0; j<=i; j++) {
				R[i][j] = R[i][j]/static_cast<double>(M);
				if(j<i)
					R[j][i] = R[i][j];
			}
		}
		// output R
		for(size_t i=0; i<n_acc; i++) {
			cout << R[i][0];
			for(size_t j=1; j<n_acc; j++) 
				cout << "\t" << R[i][j];
			cout << "\n";
		}

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

