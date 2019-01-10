#include <cmath>
#include <iostream>
#include <fstream>
#include <boost/math/distributions/chi_squared.hpp>

#include "kmers_QQ_plot_statistics.h"

using namespace std;
using boost::math::chi_squared;
using boost::math::cdf;

KmersQQPlotStatistics::KmersQQPlotStatistics():
	m_scores(0) // accumelating information (might get very large)
{}

// Just add the score to the list
void KmersQQPlotStatistics::add_score(const double &score) {
	if(score>0)  // 0 is like NaN for us
		m_scores.push_back(score);
}


///
/// @brief  print the statistics needed to paint the QQ-plot
/// @param  1. filename for output
//			2. max_points - number of sampled points from the list
//			3. actual_tests - as presence absence patterns might repeat, the actual space of test might be smaller
//
/// @return (change the order of the scores)
///
// As this class is likley to contain 10^7-10^9 data points, we won't be able to plot all of them in a QQ-plot
// Thus we will plot only a subset of the points (defined in max_points) to reach a density resonable for 
// painting on a single plot
//

void KmersQQPlotStatistics::print_stats_to_file_norm_unique_tests(const std::string &filename, const double &gamma,
		const size_t &n_individuals, const size_t actual_tests, size_t max_points) {

	float factor = 1 / (static_cast<float>(gamma) * static_cast<float>(n_individuals)); // Factor from real statistics 

	ofstream f_out(filename);
	sort(m_scores.begin(), m_scores.end());

	chi_squared dist(1); // the distribution object
	size_t counter(0);
	int index(m_scores.size()-1), step;

	// This is a small heuristic to for sampling of the scores
	// the idea is to plot all the most significant values, thus plotting the highest N/2 points, where N
	// is the total number of points we want. Then we will sample given the density of chi-square(df 1) is linear 
	// in log2 with slope y=-x. This method of sampling shold be optimized for ploting the qq-plot.
	//
	// Finding alpha is done according to this assumption & a heuristic I thought reasonable and seem work in my tests.	
	double alpha = 1;
	for(size_t i=0; i<10; i++)	
		alpha = log2((m_scores.size() - max_points/2)*(pow(2,alpha)-1)) / static_cast<double>(max_points/2);
	f_out << "chi square statistics\t-log10(p-value observed)\t-log10(p-value) expected\n";
	while((counter < max_points) && (index>0)) {
		float cur_score = m_scores[index] * factor; // as factor > 0, it is only important in the reporting
		f_out << cur_score << "\t" << -log10((1-cdf(dist, cur_score))) << "\t" << 
			-log10(norm_uniform(m_scores.size() - index, m_scores.size(), actual_tests)) << endl;

		if((2*counter)<max_points) {
			step = 1; 
		} else {
			step = pow(2,alpha * static_cast<double>(counter - max_points/2));
		}
		counter++;
		index-=step;
	}
}

// i - index / N - number of tests / K - unique number of tests (space size)
double KmersQQPlotStatistics::norm_uniform(const size_t &i, const size_t &N, const size_t &K) const {
	double i_norm = static_cast<double>((i * K+N-1) / N);
	return i_norm/static_cast<double>(K+1);
}

void KmersQQPlotStatistics::print_stats_to_file(const std::string &filename, const double &gamma,
		const size_t &n_individuals, size_t max_points) {
	print_stats_to_file_norm_unique_tests(filename, gamma, n_individuals,  m_scores.size(),
			max_points);
}
