/**
 *       @file  snps_multiple_databases.cpp
 *      @brief  Implementation of associating SNPs with phenotypes using GRAMMAR-Gamma
 *      approximation.
 *
 * Detailed description starts here.
 *
 *     @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
 *
 *   @internal
 *     Created  11/19/18
 *    Revision  $Id: doxygen.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
 *    Compiler  gcc/g++
 *     Company  Max Planck Institute for Developmental Biology Dep 6
 *   Copyright  Copyright (c) 2018, Yoav Voichek
 *
 * This source code is released for free distribution under the terms of the
 * GNU General Public License as published by the Free Software Foundation.
 * =====================================================================================
 */

#include "snps_multiple_databases.h"
#include "best_associations_heap.h"

#include <smmintrin.h>  /* SSE 4.1 */
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

using namespace std;

/* Defines */
#define ALIGNTO(n) __attribute__((aligned(n)))

/**
 *  * SSE 4.1 implementation.
 *   */
#pragma GCC push_options //This gave a small 7-8% performence improvment
#pragma GCC optimize ("unroll-loops")
inline double dot_product_SSE4(
		const std::vector<float> &S,
		std::vector<uint64_t>::const_iterator bv,
		const std::size_t &n_words
		) {
		// Variables for the scoring with SSE4
		__m128   mask ALIGNTO(16);
		__m128   f		   ALIGNTO(16);
		__m128   zblended  ALIGNTO(16);
		__m128   sums      ALIGNTO(16) = _mm_setzero_ps();
		float    sumsf[4]  ALIGNTO(16);
		
		std::size_t j=0;
		for(uint64_t hashmap_i=0; hashmap_i<n_words; hashmap_i+=2, j+=128) { // two words at a time
			mask = _mm_load_ps((const float*)&bv[hashmap_i]);
			for(size_t i=0;i<128;i+=4) {
				f = _mm_load_ps(&S[j+i]);
				zblended  = _mm_setzero_ps();					// Intialize zeros
				zblended  = _mm_blendv_ps(zblended,f, mask);	// Choose according to bits
				sums      = _mm_add_ps(sums, zblended);			// Sum chosen floats
				mask = _mm_castsi128_ps(_mm_slli_epi32(_mm_castps_si128(mask), 1));
			}
		}
		_mm_store_ps(sumsf, sums);
		return (sumsf[0] + sumsf[1] + sumsf[2] + sumsf[3]);
}
#pragma GCC pop_options




MultipleSNPsDataBases::MultipleSNPsDataBases(const string &base_name_bedbim,
		const vector<string> &samples_to_use): 
	m_base_name(base_name_bedbim),
	m_samples_names(samples_to_use),
	m_accessions_bedbim_file(),
	m_uint64_words(2*((m_samples_names.size()+(2*WLEN)-1)/(2*WLEN))), // need to multiple by 128
	m_n_snps(),
	m_n_bytes_per_snp(),
	m_presence_absence(),
	m_missing(),
	m_hetrozygous(),
	m_presence_absence_popcnt(),
	m_presence_absence_total(),
	m_S_gi_2()
{
	vector<string> all_samples = get_names_from_fam_file(m_base_name + ".fam");
	auto map_snps_table_to_current_samples = 
		create_map_from_all_samples(all_samples, m_samples_names);

	// load bed file - initial parameters
	ifstream bed_file(base_name_bedbim + ".bed", ios::binary | ios::ate); // ios::ate - put pointer in the end of file
	size_t bed_file_size = bed_file.tellg(); // file size in bytes
	bed_file.seekg(0, ios::beg);//go back to begining	
	if(bed_file_size < 3)
		throw std::logic_error("Bed file is too small");
	size_t n_samples_bedbim = all_samples.size();
	size_t n_samples = m_samples_names.size();

	m_n_bytes_per_snp = (4+n_samples_bedbim-1) / 4;
	m_n_snps = (bed_file_size-3) / m_n_bytes_per_snp;
	if(bed_file_size != ((m_n_snps * m_n_bytes_per_snp)+3)) 
		throw std::logic_error("Ilegal size of bed file");
	cerr << base_name_bedbim << "\t(snps,samples) = " << m_n_snps << ", " << n_samples_bedbim << endl;
	// read prefix
	bed_file.ignore(3);
	// build containers
	m_presence_absence.resize(m_uint64_words * m_n_snps, 0);
	m_missing.resize(m_uint64_words * m_n_snps, 0);
	m_hetrozygous.resize(m_uint64_words * m_n_snps, 0);
	
	m_presence_absence_popcnt.resize(m_n_snps, 0);
	m_presence_absence_total.resize(m_n_snps, 0);
	m_S_gi_2.resize(m_n_snps, 0);

	// read bed file into memory
	// 00 (a/a) -> 0 (+1 count | +0 popcnt)
	// 01 (NaN) -> 0 (+0 count | +0 popcnt)
	// 10 (A/a) -> 0 (+1/2 count | +1 popcnt)
	// 11 (A/A) -> 0 (+1 count | +1 popcnt)
	double _dubit_to_popcnt[] = {0,0,0.5,1};
	uint64_t _dubit_to_bit[] = {0,0,0,1};
	uint64_t _dubit_to_total[] = {1,0,1,1};
	uint64_t _dubit_to_hetrozygous[] = {0,0,1,0};

	vector<unsigned char> buffer(m_n_bytes_per_snp); // reading buffer
	size_t offset = 0;
	for(size_t i=0; i<m_n_snps; ++i) {
		bed_file.read(reinterpret_cast<char *>(buffer.data()), sizeof(char)*buffer.size());
		// Organize them into place
		double cur_popcnt = 0;
		double cur_S_gi_2 = 0;
		uint64_t cur_total = 0;
		for(size_t si=0; si<n_samples; ++si) {
			uint64_t dubit = (buffer[get<0>(map_snps_table_to_current_samples)[si]] >>
					get<1>(map_snps_table_to_current_samples)[si]) & 0x03;

			cur_popcnt += _dubit_to_popcnt[dubit];
			cur_S_gi_2 += _dubit_to_popcnt[dubit]*_dubit_to_popcnt[dubit];
			cur_total += _dubit_to_total[dubit];

			m_presence_absence[offset + (si >> 6)] ^= _dubit_to_bit[dubit] << (si & 0x3f);
			m_missing[offset + (si >> 6)] ^= _dubit_to_total[dubit] << (si & 0x3f);
			m_hetrozygous[offset + (si >> 6)] ^= _dubit_to_hetrozygous[dubit] << (si & 0x3f);
		}
		m_presence_absence_popcnt[i] = cur_popcnt;
		m_S_gi_2[i] = cur_S_gi_2;
		m_presence_absence_total[i] =  static_cast<double>(cur_total);
		offset += m_uint64_words;
	}
	// close bed file
	bed_file.close();	
}

///
/// @brief  caluclate the approximated score for associations
/// @param  phenotypes (size multiply of 128), snp index, and minor allele count (MAC)
/// @return score
///
double MultipleSNPsDataBases::calculate_grammmar_approx_association(
		const vector<float> &phenotypes,const size_t &index, const double &mac) const {
	if((mac>m_presence_absence_popcnt[index]) || (mac>(m_presence_absence_total[index]-m_presence_absence_popcnt[index]))) 
		return 0;
	size_t index_in_vector = index*m_uint64_words;
	double yigi = dot_product_SSE4(phenotypes, m_presence_absence.begin() + index_in_vector, m_uint64_words)+
		      dot_product_SSE4(phenotypes, m_hetrozygous.begin() + index_in_vector, m_uint64_words) * 0.5; // S(gi*ri)
	double score_sum = dot_product_SSE4(phenotypes, m_missing.begin() + index_in_vector, m_uint64_words); //S(vi*ri)
	double N = m_presence_absence_total[index]; // sum !is.na
	double S_gi = m_presence_absence_popcnt[index]; // S(gi) - N1
	double S_gi_2 = m_S_gi_2[index]; // S(gi*gi)

	double r =  N*yigi - S_gi*score_sum;
	r = r*r;
	return r / (N*(N*S_gi_2 -  S_gi*S_gi));
}


///
/// @brief  Reads a fam file and return the name of the samples
/// @param  
/// @return 
///
vector<string> MultipleSNPsDataBases::get_names_from_fam_file(const string &fam_fn) const {
	vector<string> sample_names;

	std::ifstream 	fin(fam_fn);
	vector<string>  line_tokens;
	string          line, cell;
	while(getline(fin, line)) { // Read the file line by line
		stringstream          lineStream(line);
		while(std::getline(lineStream, cell, ' '))
			line_tokens.push_back(cell);
		sample_names.push_back(line_tokens[0]); // second column
		line_tokens.resize(0);
	}
	return sample_names; 
}

///
/// @brief  Create a map for the 2-bits in the bits blocks creating the bed to 1 bit in our presence/absence
//	table for the subset of samples.
/// @param  	full_list - the samples from the fam file
//		sub_list -  the samples we want to keep
/// @return a map between the bed file location and our sample order
///
tuple<vector<size_t>, vector<size_t> > MultipleSNPsDataBases::create_map_from_all_samples(
		const vector<string> &full_list, const vector<string> &sub_list) const {
	tuple<vector<size_t>, vector<size_t> > res = make_tuple(vector<size_t>(sub_list.size()),
			vector<size_t>(sub_list.size()));
	for(size_t i_sub=0;i_sub<sub_list.size();++i_sub) {
		size_t i_full = find(full_list.begin(), full_list.end(), sub_list[i_sub])-
			full_list.begin();
		if(i_full == full_list.size()) 
			throw std::logic_error("All accessions should be in fam file: " + sub_list[i_sub]);
		get<0>(res)[i_sub] = (i_full / 4); // byte index
		get<1>(res)[i_sub] = (i_full % 4)*2; // dubit index in byte
	}
	return res;
}


///
/// @brief  calculate score of associations of all snps with phenotype and return the indiced of best snps
/// @param  	phenotypes values - same order as the samples name for the constructor
//		number of best snps to return
/// @return indices of snps with highest associations
///
vector<size_t> MultipleSNPsDataBases::get_most_associated_snps(vector<float> phenotypes, const size_t &number_of_best_associations, 
		const double &min_minor_allele_count) const {
	phenotypes.resize(m_uint64_words*sizeof(uint64_t)*8, 0);
	permute_scores(phenotypes);
	BestAssociationsHeap best_associations(number_of_best_associations);
	for(size_t snp_i=0; snp_i < m_n_snps; snp_i++) {
		double cur_score = calculate_grammmar_approx_association(phenotypes, snp_i, min_minor_allele_count);
		best_associations.add_association(0,cur_score, snp_i);
		//cout << snp_i << "\t" << cur_score << endl;
	}	
	return best_associations.get_rows_sorted_indices();
}



///
/// @brief  Filter snps according to their indices from the bed & bim files
/// @param  	base file names for creation of the bed & bim files
//		list of snp indices to output - NOTE: this list should be sortted! 
/// @return 
///
void MultipleSNPsDataBases::output_plink_bed_file(const vector<string> &files_base_names,
		vector<vector<size_t> > SNPs_indices) const {
	// Create output files
	vector<BedBimFilesHandle> output_files_handles;
	for(size_t i=0; i<files_base_names.size(); i++) {
		output_files_handles.emplace_back(files_base_names[i]);
	}
	// open original bed/bim files
	ifstream bed_file(m_base_name + ".bed", ios::binary);
	ifstream bim_file(m_base_name + ".bim");
	bed_file.ignore(3);
	// define buffers
	vector<char> bed_buffer(m_n_bytes_per_snp);
	string bim_buffer;
	// Last index checked for for each index list
	vector<size_t> last_index(SNPs_indices.size(), 0);
	for(size_t i=0; i<m_n_snps; i++) {
		// read snp info from file
		//bed_file.read(reinterpret_cast<char *>(bed_buffer.data()), sizeof(char)*bed_buffer.size());
		bed_file.read(bed_buffer.data(), sizeof(char)*bed_buffer.size());
		getline(bim_file, bim_buffer);
		for(size_t list_index=0; list_index<SNPs_indices.size(); list_index++) {// going over phenotypes
			if((last_index[list_index] < SNPs_indices[list_index].size()) && 
					(i == SNPs_indices[list_index][last_index[list_index]])) {

				output_files_handles[list_index].f_bim << bim_buffer << endl;
				output_files_handles[list_index].f_bed.write(bed_buffer.data(), bed_buffer.size());

				last_index[list_index]++; // advancing the index we are looking for in the list

			}
		}
	}
	// close original bed/bim files	
	bed_file.close();
	bim_file.close();

	// Close output files
	for(size_t i=0; i<output_files_handles.size(); i++)
		output_files_handles[i].close();
}

