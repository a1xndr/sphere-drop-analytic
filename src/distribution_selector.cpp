/*
 * =====================================================================================
 *
 *       Filename:  distribution_selector.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/17/2016 10:44:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <distribution_selector.hpp>

int DistributionSelector::set(dist_type dist, double a)
{
    switch(dist)
    {
	case exponential:
		d_exponential = std::exponential_distribution<double>(a);
		break;
	case chi_squared:
		d_chi_squared = std::chi_squared_distribution<double>(a);
		break;
	default:
		return -1;
    }
    type = dist;
}
int DistributionSelector::set(dist_type dist, double a , double b)
{
    switch(dist)
    {
	case normal:
		d_normal = std::normal_distribution<double> (a, b);
		break;
	case gamma:
		d_gamma = std::gamma_distribution<double>(a, b);
		break;
	case uniform:
		d_uniform = std::uniform_real_distribution<double>(a, b);
		break;
	default:
		return -1;
    }
    type = dist;
}
double DistributionSelector::gen()
{
    switch(type)
    {
	case normal:
		return d_normal(generator);
		break;
	case gamma:
		return d_gamma(generator);
		break;
	case uniform:
		return d_uniform(generator);
		break;
	case exponential:
		return d_exponential(generator);
		break;
	case chi_squared:
		return d_chi_squared(generator);
		break;
	default:
		return -1;
    }
}

