///
///      @file  kmers_QQ_plot_statistics.h
///     @brief  Class definition to hold statistics for QQ-plot
///
/// Detailed description starts here.
/// Note: The fact that we are saving all the scores is highly unneccesary, but for now this 
/// will do.
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  12/09/18
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
		KmersQQPlotStatistics();
		KmersQQPlotStatistics(const KmersQQPlotStatistics &s) = delete;
		KmersQQPlotStatistics& operator=(const KmersQQPlotStatistics) = delete;

		void add_score(const double &score);
		void print_stats_to_file(const std::string &filename, const double &gamma,
			   const std::size_t &n_individuals, std::size_t max_points = 10000);
		std::size_t total_insertions() const {return m_scores.size();}

	private:
		size_t m_n_insertions;
		std::vector<float> m_scores;
};	




#endif
