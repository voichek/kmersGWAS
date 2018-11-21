///
///      @file  kmer_general.cpp
///     @brief  This file contains general function / parameters needed for all the package
///
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

int get_mem_used_by_process_parseLine(char* line){
	// This assumes that a digit will be found and the line ends in " Kb".
	int i = strlen(line);
	const char* p = line;
	while (*p <'0' || *p > '9') p++;
	line[i-3] = '\0';
	i = atoi(p);
	return i;
}

double get_mem_used_by_process(){ //Note: this value is in MB!
	FILE* file = fopen("/proc/self/status", "r");
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != NULL){
		if (strncmp(line, "VmSize:", 7) == 0){
			result = get_mem_used_by_process_parseLine(line);
			break;
		}
	}
	fclose(file);
	return (double)result/1024.;
}


/**
 * @brief   Read a list of k-mer DBS (DB name, DB path)
 * @param   filename path
 * @return  list of DB name and path
 */
vector<KMCDataBaseHandle> read_accession_db_list(string filename) {
	ifstream fin(filename);
	vector<KMCDataBaseHandle> res;
	KMCDataBaseHandle db_info;

	while(fin >> db_info.dir_path) {
		fin >> db_info.name;
		res.push_back(db_info);
	}
	fin.close();
	return res;
}

///
/// @brief  filter from a kmer vector only the kmers part of a given set
/// @param  1. vector of uint64 representing k-mers and a set of kmers (hash set)
/// @return 
///
void filter_kmers_to_set(std::vector<uint64_t> &kmers, const KmersSet &set_kmers) {
	// need to implement....
	size_t move_to = 0;
	for(size_t i=1; i<kmers.size(); i++) {
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
