///
///      @file  snps_multiple_databases.h
///     @brief Decleration of MultipleSNPsDataBases 
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  11/18/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///
#ifndef SNPS_MULTIDB_H
#define SNPS_MULTIDB_H
#include <stdint.h>

#include <string>
#include <tuple>
#include <vector>

class MultipleSNPsDataBases {
	public:
		MultipleSNPsDataBases(const std::string &bedbin_base_fn,
				const std::vector<std::string> &samples_to_use);
		MultipleSNPsDataBases() = delete;
		MultipleSNPsDataBases(const MultipleSNPsDataBases&) = delete;
		MultipleSNPsDataBases& operator=(const MultipleSNPsDataBases&) = delete;

		~MultipleSNPsDataBases() {}; //Dtor
		// Get the indices of the best N SNPs 
		std::vector<std::size_t> get_most_associated_snps(std::vector<float> phenotypes,
				const std::size_t &number_of_best_associations,
				const double &min_minor_allele_count) const;
		// Open files for a few phenotypes together and output the given indices, respectivley
		void output_plink_bed_file(const std::vector<std::string> &files_base_names,
				std::vector<std::vector<std::size_t> > SNPs_indices) const; 
	
	private:
		std::string m_base_name;			// prefix to bedbimfam file names
		std::vector<std::string> m_samples_names;      	// Names of samples from the fam file		
		std::size_t m_accessions_bedbim_file;           // Number of samples in bedbim file

		std::size_t m_uint64_words;             	// Number of words for "row" in memory table
		std::size_t m_n_snps;				// Number of snps in file
		std::size_t m_n_bytes_per_snp;			// bytes of information in bed file per snp
		std::vector<uint64_t> m_presence_absence;       // presence/absence table
		std::vector<uint64_t> m_missing;       		// missing/non-missing table
		std::vector<uint64_t> m_hetrozygous;   		// Hetrozygous table (1 - hetrozgous | 0 - otherwise)

		std::vector<double> m_presence_absence_popcnt;	// Number of presence (hetrozygous or homozygous)
		std::vector<double> m_presence_absence_total;	// Number of presence&absence (no missing)
		std::vector<double> m_S_gi_2;			// Sum genotype square
		std::vector<std::string> get_names_from_fam_file(const std::string &fam_fn) const;
		std::tuple<std::vector<std::size_t>, std::vector<std::size_t> >  
			create_map_from_all_samples(
					const std::vector<std::string> &full_list,
					const std::vector<std::string> &sub_list) const;
		double calculate_grammmar_approx_association(const std::vector<float> &phenotypes, 
				const std::size_t &index, const double &mac) const;
};

#endif


