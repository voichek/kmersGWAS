///
///      @file  emma_kinship.cpp
///     @brief  implementation of kinship calculation from EMMA in C++
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  12/01/19
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2019, Yoav Voichek
///
///This source code is released for free distribution under the terms of the
///GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include<fstream>
#include<iostream>
#include<string>
#include<vector>

using namespace std;


void is_not_true(const bool &condition, const string &err_msg) {
	if(!condition) 
		throw std::runtime_error( "error:\t" + err_msg);
	
}

size_t count_samples_in_fam_file(const string &file_path) {
	std::ifstream 	fin(file_path);
	is_not_true(fin.is_open(),"couldn't open fam file");
	
	string          line;
	size_t counter = 0;
	while(getline(fin, line))  // Read the file line by line
		counter++;
	
	return counter; 
}
void update_K(vector<vector<double> > &K, const vector<double> &snp_calls) {
		for(size_t kr=1; kr<snp_calls.size(); ++kr) {
			for(size_t kc=0; kc<kr; ++kc) {
				K[kr][kc] += snp_calls[kr]*snp_calls[kc] + 
					(1-snp_calls[kr])*(1-snp_calls[kc]);
			}
		}
}

void plot_kinship(vector<vector<double> > &K, size_t size_plot=0xFFFFFFFF) { // assuming matrix in NxN
	if(size_plot == 0xFFFFFFFF)
		size_plot = K.size();

	for(size_t i=0; i<size_plot; ++i){
		for(size_t j=0; j<size_plot; j++) {
			if(j>0)
				cout << "\t";
			cout << K.at(i).at(j);
		}
		cout << "\n";
	}
}
vector<vector<double> > emma_kinship(const string &base_bedbim) {
	// load bed file - initial parameters
	ifstream bed_file(base_bedbim + ".bed", ios::binary | ios::ate); // ios::ate - put pointer in the end of file
	is_not_true(bed_file.is_open(), "couldn't open bed file");

	size_t bed_file_size = bed_file.tellg(); // file size in bytes
	bed_file.seekg(0, ios::beg);//go back to begining	
	is_not_true(bed_file_size >= 3,"Bed file is too small");
	
	size_t n_samples = count_samples_in_fam_file(base_bedbim + ".fam");

	size_t n_bytes_per_snp = (4+n_samples-1) / 4;
	size_t n_snps = (bed_file_size-3) / n_bytes_per_snp;
	is_not_true(bed_file_size == ((n_snps * n_bytes_per_snp)+3),"Ilegal size of bed file");

	cerr << base_bedbim << "\t(snps,samples) = " << n_snps << ", " << n_samples << endl;
	// read prefix
	bed_file.ignore(3);
	
	vector<vector<double> > K(n_samples, vector<double>(n_samples, 0));
	size_t n_snps_not_nan(0);
	for(size_t i=0; i<n_samples; i++)
		K.at(i).at(i) = 1;

	// read bed file into memory
	// 00 (a/a) -> 0 (+1 count | +0 popcnt)
	// 01 (NaN) -> 0 (+0 count | +0 popcnt)
	// 10 (A/a) -> 0 (+0 count | +0 popcnt) - we treat at hetrozygous as nan
	// 11 (A/A) -> 0 (+1 count | +1 popcnt)
	double _dubit_to_popcnt[] = {0,0,0,1};
	double _dubit_to_total[] = {1,0,1,1};
	uint64_t _dubit_is_found[] = {1,0,1,1};
	uint64_t _dubit_is_hetro[] = {0,0,1,0};
	vector<unsigned char> buffer(n_bytes_per_snp); // reading buffer
	vector<double> snp_calls(n_samples);
	vector<uint64_t> snp_not_miss(n_samples), snp_hetro(n_samples);

	for(size_t i=0; i<n_snps; ++i) {
		if((i%100000)==0) {cerr << "."; cerr.flush();}
		if((i%1000000)==0) {cerr << "M"; cerr.flush();}

		bed_file.read(reinterpret_cast<char *>(buffer.data()), sizeof(char)*buffer.size());
		// Organize them into place
		double n_var_allele = 0;
		double n_total = 0;
		uint64_t n_hetroz = 0;
		for(size_t si=0; si<n_samples; ++si) {
			uint64_t dubit = (buffer[si>>2] >> ((si%4)*2) ) & 0x03;
			n_var_allele += _dubit_to_popcnt[dubit];
			n_total += _dubit_to_total[dubit];
			snp_calls[si] = _dubit_to_popcnt[dubit];
			snp_not_miss[si] = _dubit_is_found[dubit];
			snp_hetro[si] = _dubit_is_hetro[dubit];
			n_hetroz += _dubit_is_hetro[dubit];
		}
		if(n_total>0) {
			n_snps_not_nan++;
			double maf = n_var_allele / n_total; 
			for(size_t si=0; si<n_samples; ++si) {
				if(snp_not_miss[si] == 0)
					snp_calls[si] = maf;
			}
			update_K(K, snp_calls);
			// Now change hetrozygous to 1
			n_var_allele += (double)n_hetroz;
			maf = n_var_allele / n_total;
			for(size_t si=0; si<n_samples; ++si) {
				if(snp_not_miss[si] == 0)
					snp_calls[si] = maf;
				if(snp_hetro[si] == 1) 
					snp_calls[si]=1;
			}
			update_K(K, snp_calls);
		}
	}
	cerr << endl;
	// close bed file
	bed_file.close();	
	for(size_t kr=1; kr<n_samples; ++kr) {
		for(size_t kc=0; kc<kr; ++kc) {
			K[kr][kc] /= 2.*(double)n_snps_not_nan;
			K[kc][kr] = K[kr][kc];	
		}
	}
	return K;
}

int main(int argc, char* argv[]) {
	if(argc != 2) {
		cerr << "usage: " << argv[0] << " base file name for bed/bim/fam files" << endl;
		return -1;
	}
	vector<vector<double> > K = emma_kinship(argv[1]);
	plot_kinship(K);
	return 0;
}
