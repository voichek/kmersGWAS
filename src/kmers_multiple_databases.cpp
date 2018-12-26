///
///      @file  kmer_multipleDB.cpp
///     @brief  
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  07/19/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///

#include "kmers_multiple_databases.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <cmath>
#include <nmmintrin.h> //_mm_popcnt_u64 
#include <smmintrin.h>  /* SSE 4.1 */
#include <algorithm>
#include <bitset>
#include <numeric>
#include <sstream>

using namespace std;

///
/// @brief  Ctor of MultipleKmersDataBases - 
/// @param 	DB_paths:	paths in the filesystem to directories containing the files of k-mers
//		db_names:	a list of names of all the DBs to use (same order as DB_paths)
//		sorted_kmer_fn:	filename inside each subdirectory containing the sorted k-mer list
/// @return 
///
MultipleKmersDataBases::MultipleKmersDataBases(
		const string &merge_db_file,
		const vector<string> &db_names,
		const vector<string> &db_to_use, 	
		const uint32 &kmer_len):
	m_db_names_db_file(db_names), // table in file
	m_db_names_table(db_to_use), // table in memory
	m_accessions_db_file(db_names.size()), 
	m_accessions(db_to_use.size()),
	m_kmer_table_file(merge_db_file, ios::binary | ios::ate), // ios::ate - put pointer in the end of file
	m_hash_words_db_file((m_accessions_db_file+WLEN-1)/WLEN),
	// Due to the scoring procdure, the array should be a multiplication of 128 (as two words are proccessed
	// together in the dot product calculation)
	m_hash_words(2*((m_accessions+(2*WLEN)-1)/(2*WLEN))),
	m_kmers(),
	m_kmers_table(),
	m_kmers_popcnt(),
	m_kmer_len(kmer_len),
	m_left_in_file(),
	m_kmer_number(),
	m_kmer_loaded(0),
	m_row_offset(0),
	m_map_word_index(),
	m_map_bit_index(),
	m_map_mask(),
	m_verbose(true)
{
	if(m_kmer_table_file.is_open()) {
		m_left_in_file = m_kmer_table_file.tellg(); // get file size
		if(m_left_in_file <= (4 + 8 + 4)) 
			throw std::logic_error("Kmer table size is too small");
		m_kmer_table_file.seekg (0, ios::beg); //go back to begining

		uint32 prefix, file_kmer_len;
		uint64_t file_accession_number;
		m_kmer_table_file.read(reinterpret_cast<char *>(&prefix), sizeof(prefix));
		m_kmer_table_file.read(reinterpret_cast<char *>(&file_accession_number), sizeof(file_accession_number));
		m_kmer_table_file.read(reinterpret_cast<char *>(&file_kmer_len), sizeof(file_kmer_len));
		m_left_in_file -= (sizeof(prefix) + sizeof(file_accession_number) + sizeof(file_kmer_len));

		if(prefix!=0xDDCCBBAA)
			throw std::logic_error("Incorrect prefix");
		if(file_accession_number != m_accessions_db_file)
			throw std::logic_error("Number of accession in file not as defined in class");
		if(file_kmer_len != m_kmer_len)
			throw std::logic_error("Kmer length not as defined in class");

		// From the size of the file we can calculate the number of k_mers
		size_t size_per_kmer = sizeof(uint64_t) * (1 + m_hash_words_db_file);
		if((m_left_in_file % size_per_kmer) != 0)
			throw std::logic_error("size of file not valid");
		m_kmer_number = m_left_in_file / size_per_kmer;
		cerr << "We have " << m_kmer_number << endl;
		create_map_from_all_DBs();
	} else {
		throw std::logic_error("Couldn't open kmer table file: " + merge_db_file);
	}
}


/**
 * @brief   Load k-mers from table to memory
 * The function will also squeeze to keep only k-mers which are part of the current DB set		
 * @param   size of memory batch and a set of k-mers to intersect with
 * @return return false if table file is allready finished 
 */
bool MultipleKmersDataBases::load_kmers(const uint64_t &batch_size, const KmersSet &set_kmers_to_use, const size_t &mac) {
	m_row_offset = m_kmer_loaded;
	clear();
	if(m_left_in_file == 0) {
		return false;
		m_kmer_table_file.close();
	} else { // file not empty
		vector<uint64_t> reading_buffer(m_hash_words_db_file+1);
		while((m_kmers.size()<batch_size) && (m_left_in_file>0)) {
			m_kmer_table_file.read(reinterpret_cast<char *>(reading_buffer.data()), 
					sizeof(uint64_t)*reading_buffer.size());
			//update counters
			m_left_in_file -= (reading_buffer.size() * sizeof(uint64_t));
			m_kmer_loaded++;
			if((set_kmers_to_use.size() == 0) || (lookup_x(set_kmers_to_use,reading_buffer[0]))) {
				uint64_t kmer_popcnt_orig = calculate_unsqueezed_popcnt(reading_buffer);
				if((kmer_popcnt_orig>=mac) && (kmer_popcnt_orig <= (m_accessions-mac))) {
					m_kmers.push_back(reading_buffer[0]); // save k-mer representation

					// add info to table in memory (= squeeze)
					size_t table_offset = m_kmers_table.size();
					m_kmers_table.resize(table_offset+m_hash_words, 0);
					for(size_t col_index=0; col_index<m_accessions; col_index++) {
						uint64_t new_bit = (reading_buffer[m_map_word_index[col_index]+1] 
								>> m_map_bit_index[col_index]) & 1ull;

						uint64_t hashmap_i = (col_index>>6);
						uint64_t  bit_i = col_index&63;
						m_kmers_table[table_offset + hashmap_i] |= (new_bit << bit_i);
					}
					uint64_t kmer_popcnt = 0;
					for(size_t hashmap_i=0; hashmap_i<m_hash_words; hashmap_i++)
						kmer_popcnt += _mm_popcnt_u64(m_kmers_table[table_offset+hashmap_i]);
					if(kmer_popcnt != kmer_popcnt_orig) {
						cerr << "error!!! should be equal popcnt!" << endl;
						throw std::logic_error("error");
					}
					m_kmers_popcnt.push_back((double)kmer_popcnt);
				}
			}
		}
	}
	return true;
}


uint64_t MultipleKmersDataBases::calculate_unsqueezed_popcnt(const vector<uint64_t> &table_row) const {
	uint64_t res(0);
	for(size_t i=1; i<table_row.size(); i++) // i=0 hold k-mer representation
		res += _mm_popcnt_u64(table_row[i] & m_map_mask[i-1]);
	return res;
}

///
/// @brief  output all the kmers and the presence absence info found in class
//			output to stdout (can change to output to a given stream)
/// @param  
/// @return 
///
void MultipleKmersDataBases::output_kmers_textual() const { // output all k-mers found in the hash to stdout
	// last word length in bits
	for(size_t kmer_i=0; kmer_i<m_kmers.size(); kmer_i++) {
		cout << bits2kmer31(m_kmers[kmer_i], m_kmer_len) << "\t";
		size_t container_i = kmer_i*m_hash_words;
		for(size_t i=container_i; i<(container_i+m_hash_words-1); i++) 
			cout << bitset<WLEN>(reverseOne(m_kmers_table[i])).to_string();
		cout << bitset<WLEN>(reverseOne(m_kmers_table[container_i+m_hash_words-1])).to_string() << endl;
	}
}

///
/// @brief  output all k-mers and presence/absence information in a format readable by Plink
//			that is "bed" format. We will also create a meta file of the k-mers in the bed file.
//
//			Format taken from the plink website filetypes section:
//			https://www.cog-genomics.org/plink/1.9/formats#bed
//			"The first three bytes should be 0x6c, 0x1b, and 0x01 in that order. (There are old versions of 
//			the .bed format which start with a different "magic number"; PLINK 1.9 recognizes them, but will 
//			convert sample-major files to the current variant-major format on sight. See the bottom of the 
//			original .bed definition page for details; that page also contains a more verbose version of the 
//			discussion below.)
//
//			The rest of the file is a sequence of V blocks of N/4 (rounded up) bytes each, where V is the 
//			number of variants and N is the number of samples. The first block corresponds to the first marker 
//			in the .bim file, etc.
//
//			The low-order two bits of a block's first byte store the first sample's genotype code. ("First 
//			sample" here means the first sample listed in the accompanying .fam file.) The next two bits 
//			store the second sample's genotype code, and so on for the 3rd and 4th samples. The second byte 
//			stores genotype codes for the 5th-8th samples, the third byte stores codes for the 9th-12th, etc.
//
//			The two-bit genotype codes have the following meanings:
//
//			00    Homozygous for first allele in .bim file
//			01    Missing genotype
//			10    Heterozygous
//			11    Homozygous for second allele in .bim file
//
/// @param  base path to write the .bin & .bed file
/// @return 
/// @notes	I might need to add the logic of taking out duplicates patterns here...
void MultipleKmersDataBases::output_plink_bed_file(const string &base_name) const  {
	BedBimFilesHandle f_handle(base_name);
	output_plink_bed_file(f_handle);
	f_handle.close();
};

void MultipleKmersDataBases::output_plink_bed_file(BedBimFilesHandle &f, const KmersSet &set_kmers) const  {
	for(size_t kmer_i=0; kmer_i<m_kmers.size(); kmer_i++) {
		if((set_kmers.size() == 0) || (lookup_x(set_kmers,m_kmers[kmer_i]))) { // check k-mer in set (or empty set)
			write_PA(bits2kmer31(m_kmers[kmer_i], m_kmer_len), kmer_i, f); 
		}
	}
}

void MultipleKmersDataBases::write_PA(const string &name, const size_t &kmer_i, BedBimFilesHandle &f) const {
	f.f_bim << "0\t" << name << "\t0\t0\t0\t1\n"; 
	size_t container_i = kmer_i*m_hash_words;
	size_t acc_index = 0;
	for(size_t i=container_i; i<(container_i+m_hash_words); i++) {
		uint64_t w = m_kmers_table[i];
		for(size_t bi=0; (bi<16) && (acc_index < m_accessions); bi++) { // every word is 64 accessions (16*4) every 4 accessions is a byte
			unsigned char b = (w&1);
			w >>= 1;
			b ^= ((w&1)<<2);
			w >>= 1;
			b ^= ((w&1)<<4);
			w >>= 1;
			b ^= ((w&1)<<6);
			b |= (b<<1);
			w >>= 1;
			f.f_bed << b;
			acc_index += 4;
		}
	}

}

size_t MultipleKmersDataBases::output_plink_bed_file(BedBimFilesHandle &f, const vector<AssociationOutputInfo> &kmer_list, size_t index) const {
	for(size_t kmer_i=0; kmer_i<m_kmers.size(); kmer_i++) {
		if((index<kmer_list.size()) && 
				(get<2>(kmer_list[index])==(kmer_i+m_row_offset))) {
			write_PA(bits2kmer31(get<0>(kmer_list[index]),m_kmer_len) + "_" + 
					to_string(get<1>(kmer_list[index])),
					kmer_i, f);
			index++;	
		}
	}
	return index;
}

void MultipleKmersDataBases::output_plink_bed_file_unique_presence_absence_patterns(const string &filename, KmersSet &pa_patterns) const {
	BedBimFilesHandle f_handle(filename);
	for(size_t kmer_i=0; kmer_i<m_kmers.size(); kmer_i++) {
		uint64_t seed = hash_presence_absence_pattern(kmer_i);
		if(!lookup_x(pa_patterns, seed)) { // check if new pattern!
			pa_patterns.insert(seed);
			write_PA(bits2kmer31(m_kmers[kmer_i], m_kmer_len), kmer_i, f_handle);
		}
	}
	f_handle.close();
}


/////
///// @brief  go over all the kmers now present in the multipleDB, calculate association score and add this to
////			the given heap
///// @param  kmers_and_scores - kmers heap to save the best scores
////			scores - the measurments to associate the presence/absence with
////			names_scores - name of the DB the score is relevant to
///// @return 
/////
void MultipleKmersDataBases::add_kmers_to_heap(BestAssociationsHeap &kmers_and_scores, vector<float> scores, 
		const size_t &min_cnt) const {
	// Due to efficency consideration we pass scores not by reference and change it to be a multiplication of
	// word size (saving many not neccessery "if" in the scoring procedure)
	float sum_scores = update_scores_and_sum(scores);
	for(size_t kmer_index=0; kmer_index<m_kmers.size(); kmer_index++)
		kmers_and_scores.add_association(m_kmers[kmer_index], 
				calculate_kmer_score(kmer_index, scores, sum_scores,  min_cnt),
				m_row_offset + kmer_index);
}


// Also save statistics on the scores, for plotting a QQ plot
void MultipleKmersDataBases::add_kmers_to_heap(BestAssociationsHeap &kmers_and_scores, vector<float> scores, 
		const size_t &min_cnt, KmersQQPlotStatistics &qq_plot_stats) const {
	float sum_scores = update_scores_and_sum(scores);
	for(size_t kmer_index=0; kmer_index<m_kmers.size(); kmer_index++) {
		double current_kmer_stat = calculate_kmer_score(kmer_index, scores, sum_scores,  min_cnt);
		kmers_and_scores.add_association(m_kmers[kmer_index], current_kmer_stat, m_row_offset + kmer_index);
		qq_plot_stats.add_score(current_kmer_stat); // saving statistics on qq-plot
	}
}

float MultipleKmersDataBases::update_scores_and_sum(vector<float> &scores) const {
	scores.resize(m_hash_words*sizeof(uint64_t)*8, 0); 
	permute_scores(scores);
	float sum_scores(0);
	for(size_t i=0; i<scores.size(); i++) 
		sum_scores += scores[i];
	return sum_scores;
}
// return the indices of DB names inserted in the class DBs
void MultipleKmersDataBases::create_map_from_all_DBs() {
	m_map_word_index.resize(0);
	m_map_bit_index.resize(0);
	m_map_mask.resize(m_hash_words_db_file, 0);
	size_t index_file_table;
	for(size_t index_mem_table=0; index_mem_table < m_accessions; index_mem_table++) {
		index_file_table = find(m_db_names_db_file.begin(), m_db_names_db_file.end(), m_db_names_table[index_mem_table]) 
			- m_db_names_db_file.begin();
		if(index_file_table == m_accessions_db_file) // Couldn't find it
			throw std::logic_error("All accessions suppose to be in DB file: " + m_db_names_table[index_mem_table]);
		m_map_word_index.push_back(index_file_table / WLEN);
		m_map_bit_index.push_back(index_file_table % WLEN);
		m_map_mask[m_map_word_index.back()] |= (1ull << m_map_bit_index.back());
	}
}


///     Very efficent dot product between bit-array (uint64_t) to floats array
///	using the specific archtecture of the proccessor (SSE4). 
///	Taken from: https://stackoverflow.com/questions/16051365/fast-dot-product-of-a-bit-vector-and-a-floating-point-vector
///	There is another version using the AVX2 archtecture which is 2 times faster.
///	But part of our cluster in MPI can't run AVX2 so we will stick to generality for now.
/* Defines */
#define ALIGNTO(n) __attribute__((aligned(n)))

/**
 *  * SSE 4.1 implementation.
 *   */
#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
double MultipleKmersDataBases::calculate_kmer_score(
		const size_t kmer_index, 
		const vector<float> &scores, // Scores has to be a multiplication of wordsize (64) 
		const float score_sum,
		const uint64_t min_in_group
		) const {

	double N = (double)m_accessions;
	double N1 = m_kmers_popcnt[kmer_index];
	double N0 = N-N1;
	if((min_in_group <= N0) && (min_in_group <= N1)) {
		// Variables for the scoring with SSE4
		__m128   mask 		ALIGNTO(16);
		__m128   f	  	ALIGNTO(16);
		__m128   zblended 	ALIGNTO(16);
		__m128   sums     	ALIGNTO(16) = _mm_setzero_ps();
		float    sumsf[4] 	ALIGNTO(16);

		size_t container_i = kmer_index*m_hash_words;
		size_t j=0;
		for(uint64_t hashmap_i=0; hashmap_i<m_hash_words; hashmap_i+=2, j+=128) { // two words at a time
			mask = _mm_load_ps((const float*)&m_kmers_table[container_i+hashmap_i]);
			for(size_t i=0;i<128;i+=4) {
				f = _mm_load_ps(&scores[j+i]);
				zblended  = _mm_setzero_ps();			// Intialize zeros
				zblended  = _mm_blendv_ps(zblended,f, mask);	// Choose what to sum according to bits
				sums      = _mm_add_ps(sums, zblended);		// Sum chosen floats
				mask = _mm_castsi128_ps(_mm_slli_epi32(_mm_castps_si128(mask), 1)); //shift
			}
		}
		_mm_store_ps(sumsf, sums);
		double yigi = sumsf[0] + sumsf[1] + sumsf[2] + sumsf[3];
		double r =  N*yigi- N1*score_sum;
		r = r * r;
		return r / (N*N1 -  N1*N1);
	} else {return 0;}
}
#pragma GCC pop_options

// Calculate hash of presence absence pattern
uint64_t MultipleKmersDataBases::hash_presence_absence_pattern(const size_t &kmer_index) const {
	static Hash64 hasher;
	size_t container_i = kmer_index*m_hash_words;
	uint64_t seed(0);
	for(uint64_t hashmap_i=0; hashmap_i<m_hash_words; hashmap_i++) 
		seed ^= hasher(m_kmers_table[container_i+hashmap_i]) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
	return seed;
}

//Saving hash of presence/absence patterns
void MultipleKmersDataBases::update_presence_absence_pattern_counter(KmersSet &pa_pattern_counter) const {
	for(size_t kmer_index=0; kmer_index<m_kmers.size(); kmer_index++) 
		pa_pattern_counter.insert(hash_presence_absence_pattern(kmer_index));// update set	
}


///
/// @brief  To calculate the gamma correction factor parameters calculated on genotypes need to be used
//			this function calculate this params.
/// @param  R -	is the matrix of the accumelated calculations
//			M -	Counter of the number of genotypes updated into R
/// @return Just update the input parameters
///
void MultipleKmersDataBases::update_gamma_precalculations(vector<vector<double> > &R, size_t &M) const {
	double n  = (double)m_accessions;
	vector<double> g(m_accessions, 0);

	for(size_t kmer_index=0; kmer_index<m_kmers.size(); kmer_index++) {
		double Egm = m_kmers_popcnt[kmer_index] / n;
		double denominator_factor = sqrt(1/(n*(Egm-Egm*Egm)));

		// Calculate g's
		size_t container_i = kmer_index*m_hash_words;
		for(size_t i=0; i<m_accessions; i++) {
			uint64_t hashmap_i = (i>>6);
			uint64_t bit_i = i&63;
			g[i] = static_cast<double>((m_kmers_table[container_i+hashmap_i] >> bit_i) & 1ull)-Egm;
			g[i] *= denominator_factor;
		}

		// update R
		for(size_t i=0; i<m_accessions; i++) {
			for(size_t j=0; j<=i; j++) {
				R[i][j] += g[i]*g[j];
			}
		}
		// update M
		M++;
	}
}

void MultipleKmersDataBases::update_emma_kinshhip_calculation(std::vector<std::vector<uint64_t> > &K, 
		uint64_t &counter) const {
	vector<uint64_t> g(m_accessions, 0);
	for(size_t kmer_index=0; kmer_index<m_kmers.size(); kmer_index++) {
		// Calculate g's
		size_t container_i = kmer_index*m_hash_words;
		for(size_t i=0; i<m_accessions; i++) {
			uint64_t hashmap_i = (i>>6);
			uint64_t bit_i = i&63;
			g[i] = (m_kmers_table[container_i+hashmap_i] >> bit_i) & 1ull;
		}

		// update R
		for(size_t i=0; i<m_accessions; i++) 
			for(size_t j=0; j<i; j++) 
				K[i][j] += (1ull^g[i]^g[j]);

		// update M
		counter++;
	}
}


vector<vector<double> > read_matrix_tab_seperated(const string &filename) {
	ifstream fin(filename);
	vector<string>  line_tokens;
	string          line, cell;

	size_t line_n(0);
	size_t N(0); //matrix row/column dimension
	
	vector<vector<double> > M;
	while(getline(fin, line)) { // Read the file line by line
		stringstream          lineStream(line);
		while(std::getline(lineStream, cell, '\t'))
			line_tokens.push_back(cell);
		if(line_n == 0)	
			N = line_tokens.size();
		if(line_tokens.size() != N) 
				throw std::logic_error("File should have the same number of fields in each row - " + filename);
		M.resize(line_n+1);
		for(size_t i=0; i<N; i++) M[line_n].push_back(stof(line_tokens[i])); 
		line_n++;
		line_tokens.resize(0);
	}
	if(line_n != N)
		throw std::logic_error("File should have the same number of raws & columns - " + filename);
	return M;
}

double calc_gamma(const string &fn_inv_cov_mat, MultipleKmersDataBases &multiDB, 
		const size_t &min_count,	const size_t max_variants) {
	size_t M(0);
	vector<vector<double> > inv_cov_mat = read_matrix_tab_seperated(fn_inv_cov_mat);
	size_t n_acc = inv_cov_mat.size();
	vector<vector<double> > R(n_acc, vector<double>(n_acc, 0));

	size_t batch_size = 10000; // constant for now
	size_t max_batch = (max_variants+batch_size) / batch_size;
	size_t batch_index = 0;
	while(multiDB.load_kmers(batch_size, min_count) && (batch_index < max_batch)) {
		multiDB.update_gamma_precalculations(R, M);
		batch_index++;
	}

	for(size_t i=0; i<n_acc; i++) {
		for(size_t j=0; j<=i; j++) {
			R[i][j] = R[i][j]/static_cast<double>(M);
			if(j<i)
				R[j][i] = R[i][j];
		}
	}

	double gamma = 0;
	for(size_t i=0; i<n_acc; i++)
		for(size_t j=0; j<n_acc; j++)
			gamma += inv_cov_mat[i][j] * R[i][j];

	return gamma;
}
