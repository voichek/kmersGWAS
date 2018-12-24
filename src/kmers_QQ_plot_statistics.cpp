#include <cmath>
#include <iostream>
#include <fstream>
#include <boost/math/distributions/chi_squared.hpp>

#include "kmers_QQ_plot_statistics.h"

using namespace std;
using boost::math::chi_squared;
using boost::math::cdf;

KmersQQPlotStatistics::KmersQQPlotStatistics(const double &gamma, const double &N_individuals):
	m_factor(1 / (gamma * N_individuals)), // Factor from real statistics 
	m_scores(0) // accumelating information (might get very large)
{}

// Just add the score to the list
void KmersQQPlotStatistics::add_score(double score) {
	if(score>0)  // 0 is like NaN for us
		m_scores.push_back(score*m_factor);
}


///
/// @brief  print the statistics needed to paint the QQ-plot
/// @param  1. filename for output
//			2. max_points - number of sampled points from the list
/// @return (change the order of the scores)
///
// As this class is likley to contain 10^7-10^9 data points, we won't be able to plot all of them in a QQ-plot
// Thus we will plot only a subset of the points (defined in max_points) to reach a density resonable for 
// painting on a single plot
//
void KmersQQPlotStatistics::print_stats_to_file(const std::string &filename, size_t max_points) {
	ofstream f_out(filename);
	sort(m_scores.begin(), m_scores.end());

	chi_squared dist(1); // the distribution object
	size_t counter(0);
	int index(m_scores.size()-1), step;

	// This is a small heuristic to the sampling of the scores
	// the idea is to plot all the most significant values, thus plotting the highest N/2 points, where N
	// is the total number of points we want. Then we will sample given the density of chi-square(df 1) is linear 
	// in log2 with slope y=-x. This method of sampling shold be optimized for ploting the qq-plot.
	//
	// Finding alpha is done according to this assumption and a heuristic I thought reasonable and seem work in my tests.	
	double alpha = 1;
	for(size_t i=0; i<10; i++)	
		alpha = log2((m_scores.size() - max_points/2)*(pow(2,alpha)-1)) / static_cast<double>(max_points/2);
	f_out << "chi square statistics\t-log10(p-value observed)\t-log10(p-value) expected\n";
	while((counter < max_points) && (index>0)) {
		f_out << m_scores[index] << "\t" << -log10((1-cdf(dist, m_scores[index]))) << "\t" << 
			-log10(static_cast<double>(m_scores.size() - index)/static_cast<double>(m_scores.size()+1)) << endl;

		if((2*counter)<max_points) {
			step = 1; 
		} else {
			step = pow(2,alpha * static_cast<double>(counter - max_points/2));
		}
		counter++;
		index-=step;
	}
}
