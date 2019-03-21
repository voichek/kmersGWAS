///
///      @file  kmers_table_to_bed.cpp
///     @brief  Take a kmers table and convert it to bed/bim file(s)
///
/// In some case we would like to run the full kmers presence/absence information in the exact
/// LMM method of choice (e.g. GEMMA). For this purpose we need to convert the kmers table to 
/// bed/bim/fam format, which is the purpose of this code. To save space, we will filter the 
/// kmers to have the wanted minor allele frequency / count (MAF, MAC) before outputing them.
/// Moreover, we will also filter kmers which have the same exact presence/absence pattern.
/// This will be done to reduce the computational load on the next step (e.g. GEMMA).
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

#include "kmer_general.h"
#include "kmers_multiple_databases.h"


using namespace std;

int main(int argc, char* argv[]) {
	if(argc != 8) {
		cerr << "usage: " << argv[0] << " <kmers table> <phenotype file> <MAF> <MAC> <max batch size> <kmer_len> <output base name>" << endl;
		cerr << "\t\tIf there is more than one phenotype in file, take the first one to the fam" << endl;
		return -1;
	}
	
	// Load parameters
	string fn_kmers_table(argv[1]);
	string fn_phenotypes(argv[2]);
	double MAF(atof(argv[3]));
	size_t MAC(atoi(argv[4]));
	size_t max_batch_size(atoi(argv[5]));
	size_t kmer_len(atoi(argv[6]));
	string output_base(argv[7]);

	
	// load phenotypes (needed for current DBs)
	pair<vector<string>, vector<PhenotypeList>> phenotypes_file_info = load_phenotypes_file(fn_phenotypes);
	cerr << "using " << phenotypes_file_info.first[0] << endl;
	PhenotypeList pheno_info = intersect_phenotypes_to_present_DBs(phenotypes_file_info.second[0], fn_kmers_table, true);
		
	// calculate effective MAC (if MAF indicate lower number)
	size_t min_count = ceil(double(pheno_info.first.size())*MAF); // MAF of 5% - maybe should make this parameter external
	if(min_count < MAC)
		min_count = MAC;

	// create kmers-multiple dbs file
	MultipleKmersDataBases multiDB(
			fn_kmers_table,
			pheno_info.first,
			kmer_len);


	// create presence/absence pattern counter (initialize with 100M)
	KmersSet pa_patterns_counter(1000*1000*100); 
	pa_patterns_counter.set_empty_key(NULL_KEY); // need to define empty value for google dense hash table

	size_t batch_index(0);
	cerr << "loading.... " << endl;
	while(multiDB.load_kmers( max_batch_size, min_count)) { // Running over the kmers table in batches
		cerr << "Batch:\t" << batch_index+1 << endl;
		// output current batch
		string fn_bedbim = output_base + "." + std::to_string(batch_index);
		multiDB.output_plink_bed_file_unique_presence_absence_patterns(fn_bedbim, pa_patterns_counter);
		// output fam file (the same in all repeats)
		ofstream fout(fn_bedbim + ".fam");
		for(size_t i=0; i<pheno_info.first.size(); i++) 
			fout << pheno_info.first[i] << " " <<
				pheno_info.first[i] << " 0 0 0 " <<
				pheno_info.second[i] << endl;
		
		fout.close();
		cerr << "done. Total parsed k-mers:\t " << pa_patterns_counter.size() << endl;
		batch_index++; // advance batch counter
	}	

	return 0;
}
