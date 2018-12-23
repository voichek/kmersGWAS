#include "kmers_QQ_plot_statistics.h"
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
KmersQQPlotStatistics::KmersQQPlotStatistics(const double &gamma, const double &N_individuals, const double &resolution, const double &max_statistics):
	m_n_insertions(0), // counter
	m_resolution(resolution), // in which resoultion to save the numbers
	m_factor(1 / (gamma * N_individuals)), // Factor from real statistics 
	m_max_value(max_statistics), // Maximal possible value of stats to save
	m_stats(ceil(m_max_value / m_resolution)+1,0) 
{}

void KmersQQPlotStatistics::add_score(double score) {
	if(score>0) {
		score *= m_factor;
	//	cerr << score << endl;
		m_n_insertions++;
		if(score>m_max_value)
			score = m_max_value;
		size_t i = ceil(score / m_resolution);
		m_stats[i]++;
	}
}

void KmersQQPlotStatistics::print_stats_to_file(const std::string &filename) const {
	ofstream f_out(filename);
	for(size_t i=0; i<m_stats.size(); i++) {
		f_out << i << "\t" << ((i+0.5)*m_resolution) << "\t" << m_stats[i] << endl;
	}
	f_out.close();
}
