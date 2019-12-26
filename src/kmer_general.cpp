///
///      @file  kmer_general.cpp
///     @brief  This file contains general function / parameters needed for all the package
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  10/28/18
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmer_general.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <sstream>

using namespace std;

/**
 * @brief   Read a list of k-mer DBS (DB name, DB path)
 * @param   filename path
 * @return  list of DB name and path
 */
vector<AccessionPath> read_accessions_path_list(string filename) {
	ifstream fin(filename);
	vector<AccessionPath> res;
	AccessionPath db_info;

	while(fin >> db_info.path) {
		fin >> db_info.name;
		res.push_back(db_info);
	}
	fin.close();
	return res;
}

vector<string> load_kmers_talbe_column_names(const string &kmers_table_base) {
	ifstream fin(kmers_table_base + string(".names"));
	vector<string> res;
	string word;
	while(fin >> word)
		res.push_back(word);
	fin.close();
	return res;
}

///
/// @brief  filter from a kmer vector only the kmers part of a given set
/// @param  vector of uint64 representing k-mers and a set of kmers (hash set)
/// If kmers were sorted they would stay sorted
///
void filter_kmers_to_set(std::vector<uint64_t> &kmers, const KmersSet &set_kmers) {
	size_t move_to = 0;
	for(size_t i=0; i<kmers.size(); i++) {
		if(lookup_x(set_kmers,kmers[i])) {
			kmers[move_to] = kmers[i];
			move_to++;
		}
	}
	kmers.resize(move_to);
}


/**
 * @brief  Given a k-mer representation in bits convert it to string
 * @param  64-bit k-mer representation 
 * @return  string with bp representation
 */
string bits2kmer31(uint64_t w, const size_t& k) {
	const static char dict_bp[] = {'A','C','G','T'};
	const static uint64_t mask2bits = 0x0000000000000003;

	string res(k,'X');
	for(std::size_t i=0; i<k; i++) {
		res[k-1-i] = dict_bp[w & mask2bits];
		w = (w>>2);
	}
	return res;
}

uint64_t reverseOne(uint64_t x) {
	x = ((x & 0xFFFFFFFF00000000) >> 32) | ((x & 0x00000000FFFFFFFF) << 32);
	x = ((x & 0xFFFF0000FFFF0000) >> 16) | ((x & 0x0000FFFF0000FFFF) << 16);
	x = ((x & 0xFF00FF00FF00FF00) >> 8)  | ((x & 0x00FF00FF00FF00FF) << 8);
	x = ((x & 0xF0F0F0F0F0F0F0F0) >> 4)  | ((x & 0x0F0F0F0F0F0F0F0F) << 4);
	x = ((x & 0xCCCCCCCCCCCCCCCC) >> 2)  | ((x & 0x3333333333333333) << 2);
	x = ((x & 0xAAAAAAAAAAAAAAAA) >> 1)  | ((x & 0x5555555555555555) << 1);
	return x;
}

/**
 * @brief   Return the time in seconds (for measuring performence)
 */
double get_time(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + (tv.tv_usec / 1000000.0);
}


/**
 * @brief   load a file that have a list of k-mers, with or without scores 
 * @param  file path, initial size of hash set and a flag indicating if the file is with or without scores
 * @return  KmersSet (hash set) having all the k-mers (without scores)
 */
KmersSet load_kmer_raw_file(string filename, size_t set_initial_size, const bool with_scores) {
	KmersSet kmer_list_to_use(set_initial_size);
	kmer_list_to_use.set_empty_key(NULL_KEY); // need to define empty value for google dense hash table

	ifstream kmer_file(filename, std::ifstream::binary);
	if(kmer_file) { // if file could be open
		kmer_file.seekg(0, kmer_file.end); // 
		uint64_t kmers_in_file = (kmer_file.tellg()) >> 3;
		if(with_scores) // Two words for each k-mer
			kmers_in_file = (kmers_in_file >> 1);
		kmer_file.seekg(0, kmer_file.beg);
		uint64_t kmer_uint;
		double score;
		if(with_scores) {
			for(uint64_t i=0; i<(kmers_in_file); i++) {
				kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
				kmer_file.read(reinterpret_cast<char *>(&score), sizeof(score));
				kmer_list_to_use.insert(kmer_uint);
			}
		} else {
			for(uint64_t i=0; i<(kmers_in_file); i++) {
				kmer_file.read(reinterpret_cast<char *>(&kmer_uint), sizeof(kmer_uint));
				kmer_list_to_use.insert(kmer_uint);
			}
		}
		kmer_file.close();
	}
	return kmer_list_to_use;
}

// Close handles of bed & bim files
void BedBimFilesHandle::close() {
	if(f_bed.is_open())
		f_bed.close();
	if(f_bim.is_open())
		f_bim.close();
}


// Permuting the order of the phenotypes for the implementation of the scoring
void permute_scores(std::vector<float> &V) { // assume V is a multiplication of 128
	std::vector<float> R(V.size());
	size_t index =0;
	for(size_t offset=0; offset<V.size(); offset+=128) {
		for(size_t i=0; i<32; i++) {
			for(size_t j=0; j<128; j+=32) {
				R.at(index) = V.at(31-i+j+offset);//at just for checking 
				index++;
			}
		}
	}
	V.swap(R);
}

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

size_t get_index_DB(const string &name, const vector<string> &names) {
	size_t index = (~0u);
	for(size_t j=0; j<names.size(); j++) {
		if(names[j] == name) {
			if(index != (~0u))
				throw std::logic_error("Two DBs with the same name! " + name);
			index = j;
		}
	}
	return index;
}

PhenotypeList intersect_phenotypes_to_present_DBs(const PhenotypeList &pl, const string &kmers_table_base, const bool &must_be_present) {
	PhenotypeList intersect_pl;
	vector<string> accessions_names = load_kmers_talbe_column_names(kmers_table_base); 
	for(size_t i=0; i<pl.first.size(); i++) {
		size_t index = get_index_DB(pl.first[i], accessions_names);
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

uint64_t kmers_step_to_threshold(const uint64_t &step, const uint64_t &total_steps, const uint64_t &kmer_length) {
	uint64_t max_kmer = ((1ull << (kmer_length*2ull))-1ull); // = 00..00111..111
	return ((max_kmer / total_steps)+1)*step;
}

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
			default : throw std::logic_error("Ilegal kmer");
					  break;
		}
		b |= (cur_dubit << (i*2));
	}
	uint64_t bt = kmer_reverse_complement(b, k.size());  

	if (bt < b)
		return bt;
	else
		return b;
}

bool is_file_exist(const char *fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}

bool is_file_exist(const string &fileName)
{
	std::ifstream infile(fileName);
	return infile.good();
}
