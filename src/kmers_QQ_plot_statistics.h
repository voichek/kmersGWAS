///
///      @file  kmers_QQ_plot_statistics.h
///     @brief  Class definition to hold statistics for QQ-plot
///
/// Detailed description starts here.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  12/09/18
///   Revision  $Id: doxygen.cpp.templates,v 1.3 2010/07/06 09:20:12 mehner Exp $
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2018, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///
#ifndef KMERS_QQ_PLOT_STATS_H 
#define KMERS_QQ_PLOT_STATS_H 

#include <string>
#include <vector>

class KmersQQPlotStatistics {
	public:
		KmersQQPlotStatistics(const double &gamma, 
				const double &N_individuals);

		KmersQQPlotStatistics() = delete;
		KmersQQPlotStatistics(const KmersQQPlotStatistics &s) = delete;
		KmersQQPlotStatistics& operator=(const KmersQQPlotStatistics) = delete;

		void add_score(const double score);
		void print_stats_to_file(const std::string &filename, std::size_t max_points = 10000);
		std::size_t total_insertions() const {return m_scores.size();}

	private:
		size_t m_n_insertions;
		double m_factor;
		std::vector<float> m_scores;
};	




#endif
