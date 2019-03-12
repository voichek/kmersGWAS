///
///      @file  filter_kmers.cpp
///     @brief  Given a list of k-mers and a k-mers table output their presence/absence pattern
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  01/14/19
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2019, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include "kmer_general.h" // read DBs

using namespace std;

uint64_t kmer2bits(string k) {
	uint64_t b(0);
	uint64_t cur_dubit;
	for(size_t i=0; i<k.size(); i++) {
		switch(char(k[k.size()-i-1])) {
			case 'A' :cur_dubit = 0 ;
			        break;
			case 'C': cur_dubit = 1 ;
			       break;
			case 'G': cur_dubit = 2 ;
			       break;
			case 'T': cur_dubit = 3 ;
				  break;
			default: throw std::logic_error("Ilegal kmer");
				 break;
		}
		b |= (cur_dubit << (i*2));
	}
//	uint64_t bt = (~reverseOne(b)) >> (WLEN-(k.size()*2));
//	
//	if (bt < b)
//		return bt;
//	else
		return b;
}

//check all k-mers are from the same size
pair<vector<uint64_t>, uint32_t> read_and_sort_kmers(const string file_name) {
	string word;
	size_t index(0);
	uint32_t kmer_len(0);
	vector<uint64_t> kmer_list;
	ifstream fin(file_name);
	while(fin >> word) {
		if(index == 0) {
			kmer_len = word.size();
		}
		if(word.size() != kmer_len) {
			cerr << "all kmers should be of the same size: " << word << endl;
			throw std::logic_error("kmers of different size");
		}
		kmer_list.push_back(kmer2bits(word));
		//cout << word << "\t" << kmer2bits(word) << "\t" << bits2kmer31(kmer2bits(word),31) <<  endl;
		index++;
	}
	fin.close();
	// sort k-mers
//	sort(kmer_list.begin(), kmer_list.end());
	return make_pair(kmer_list, kmer_len);
}

int main(int argc, char* argv[]) {
	if(argc != 5) {
		cerr << "usage: " << argv[0] << " <file with kmers list> <kmers table> <DB file> <output file>" << endl;
		return 1;
	}

	vector<uint64_t> sorted_kmers;
	uint32_t kmer_len;
	tie(sorted_kmers, kmer_len) = read_and_sort_kmers(argv[1]);
	if(sorted_kmers.size() == 0) {
		cerr << "kmers file is empty" << endl;
		return 1;
	}
	vector<AccessionPath> DB_paths = read_accessions_path_list(argv[3]);
	uint64_t words_per_kmer = (DB_paths.size() +  WLEN - 1) / WLEN;

	ifstream table_handle(argv[2], ios::binary | ios::ate);
	if(table_handle.is_open()) {
		size_t left_in_file = table_handle.tellg();
		if(left_in_file <= (4 + 8 + 4)) {
			cerr << "table file is too small" << endl;
			return 1;
		}
		table_handle.seekg(0, ios::beg); // go to file begining

		uint32 prefix, file_kmer_len;
		uint64_t file_accession_number;
		table_handle.read(reinterpret_cast<char *>(&prefix), sizeof(prefix));
		table_handle.read(reinterpret_cast<char *>(&file_accession_number), sizeof(file_accession_number));
		table_handle.read(reinterpret_cast<char *>(&file_kmer_len), sizeof(file_kmer_len));
		left_in_file -= (sizeof(prefix) + sizeof(file_accession_number) + sizeof(file_kmer_len));

		if(prefix!=0xDDCCBBAA)
			throw std::logic_error("Incorrect prefix");
		if(file_accession_number != DB_paths.size() )
			throw std::logic_error("number of accession in file not as defined in class");
		if(file_kmer_len != kmer_len)
			throw std::logic_error("kmer length in table and in list are not the same");

		// From the size of the file we can calculate the number of k_mers
		size_t size_per_kmer = sizeof(uint64_t) * (1 + words_per_kmer);
		if((left_in_file % size_per_kmer) != 0)
			throw std::logic_error("size of file not valid");

		uint64_t kmer_number = left_in_file / size_per_kmer;
		cerr << "We have " << kmer_number << endl;

		// Open output file
		ofstream fout(argv[4]);
		if(!fout.is_open()) {
			cerr << "can't open output file " << endl;
			return 1;
		}
		fout << "kmer"; // output header
		for(size_t i=0; i<DB_paths.size(); i++)
			fout << "\t" << DB_paths[i].name;
		fout << "\n";

		// Start reading files
		uint64_t i_kl(0); // index kmers list
		uint64_t i_kt(0); // index kmers table
		vector<uint64_t> buffer(words_per_kmer+1);
		while((i_kl < sorted_kmers.size()) && (i_kt < kmer_number)) {
			table_handle.read(reinterpret_cast<char *>(buffer.data()), 
					sizeof(uint64_t)*buffer.size());
			if(buffer[0] == sorted_kmers[i_kl]) { // k-mer to find
				//cerr << sorted_kmers[i_kl] << endl;
				fout << bits2kmer31(buffer[0], kmer_len);
				//output presence/ absence	
				for(size_t col_index=0; col_index<DB_paths.size(); col_index++) {
					uint64_t new_bit = buffer[(col_index >> 6) + 1] >> (col_index&(WLEN-1))  & 1ull;
					fout << "\t" << new_bit;
				}
				fout << "\n";
				i_kl++; // next k-mer to look for
			}
			i_kt++; // next row in table
		}
		fout.close();
		table_handle.close();
	} else {
		cerr << "Can't open table file" << endl;
		return 1;
	}

	return 0;
}
