///
///      @file  best_associations_map.h
///     @brief  Decleration of k-mers heap 
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  11/15/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///


#ifndef BEST_ASSOCIATIONS_H
#define BEST_ASSOCIATIONS_H

#include "kmer_general.h"

/***********************************************************************************************************/
/**
 * @class BestAssociationsHeap
 * @brief save a priority score of kmers
 */
class BestAssociationsHeap {
	public:
		BestAssociationsHeap(std::size_t max_results);
		void add_association(const uint64_t &k, const double &score, const uint64_t &kmer_row);

		void output_to_file(const std::string &filename) const;
		void output_to_file_with_scores(const std::string &filename) const;
		// it would be nice to some how create a histogram of all the scores along the way...
		// (not only the ones we keep)
		void plot_stat() const;
		inline void empty_heap() {AssociationsPriorityQueue().swap(m_best_kmers);} // empty heap content
		KmersSet get_KmersSet() const;
		kmers_output_list get_kmers_for_output(const size_t &kmer_len) const;
		std::vector<std::size_t> get_rows_sorted_indices() const;
		std::size_t number_of_insertion() const {return cnt_kmers;}
	private:
		std::size_t m_n_res;
		AssociationsPriorityQueue m_best_kmers; // heap that will contain the scores
		size_t cnt_kmers;
		size_t cnt_pops;
		size_t cnt_push;
		double lowest_score;	
};

#endif
