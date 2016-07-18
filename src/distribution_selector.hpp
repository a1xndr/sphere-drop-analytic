/*
 * =====================================================================================
 *
 *       Filename:  distribution_selector.hpp
 *
 *    Description:  A mux for several std distribution classes
 *        Version:  1.0
 *        Created:  07/17/2016 10:11:36 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexander Oleinik
 *   Organization:  
 *
 * =====================================================================================
 */
#include <random>

class DistributionSelector
{
    public:
        enum dist_type{
            normal,
            gamma,
            uniform,
            exponential,
            chi_squared
        };
        dist_type type;
        int set(dist_type dist, double a);
        int set(dist_type dist, double a , double b);
        double gen();
    private:
        std::normal_distribution<double>        d_normal;
        std::gamma_distribution<double>         d_gamma;
        std::uniform_real_distribution<double>       d_uniform;
        std::exponential_distribution<double>   d_exponential;
        std::chi_squared_distribution<double>   d_chi_squared;
        std::default_random_engine generator;

};

