#ifndef MEAN_POP_DIST_H
#define	MEAN_POP_DIST_H

#include <mutation_accumulation/utility/distribution_traits.h> // distribution_types::pmf_type
#include <mutation_accumulation/utility/data_traits.h> // data_types::discrete_type
#include <mutation_accumulation/statistics/statistics/distribution_statistics_population.h> // monte_carlo::Distribution_Statistics_Population

/*************************************************************************/

namespace monte_carlo {

    namespace mean_detail {

        /** 
         * calculate mean given that distribution is PMF, \n
         * which implies that sample type is discrete: int, long long int, bool \n
         * assumes that PMF is defined for each integer \n
         * assumes that pmf reaches close enough to zero
         * TESTED
         */
        template<class sample_type>
        const double calculate_mean(
        const probability::SampleSpace_Probability<sample_type> &sampleSpace_probability,
        distribution_types::pmf_type) {

            const std::vector<sample_type> _sample_space = sampleSpace_probability.sample_space();
            const std::vector<double> _probability = sampleSpace_probability.probability();

            double mean = 0.0;
            for (int ii = 0; ii < _sample_space.size(); ii++)
                mean += static_cast<double> (_sample_space.at(ii)) * _probability.at(ii);

            return mean;
        }

        /** 
         * calculate mean given that distribution is PDF \n
         * <X> \approx \sum_bins bin_center * P(X in bin) \n
         * this function is yet to be implemented
         */

        /** 
         * calculate mean given that distribution is CDF \n
         * sample type is discrete or continuous\n
         */
        template<class sample_type>
        const double calculate_mean(
        const probability::SampleSpace_Probability<sample_type> &sampleSpace_probability,
        distribution_types::cdf_type) {

            typename data_types::data_traits<sample_type>::category sample_category;
            return calculate_mean_cdf(sampleSpace_probability, sample_category);
        }

        /** 
         * calculate mean given that distribution is CDF and sample type is discrete\n
         * assumes that CDF is defined for each integer \n
         * assumes that cdf reaches close enough to unity
         * TESTED
         */
        template<class sample_type>
        const double calculate_mean_cdf(
        const probability::SampleSpace_Probability<sample_type> &sampleSpace_probability,
        data_types::discrete_type) {

            const std::vector<double> _probability = sampleSpace_probability.probability();

            double mean = 0.0;
            for (int ii = 0; ii < _probability.size(); ii++)
                mean += 1.0 - _probability.at(ii);

            return mean;
        }

        /** 
         * calculate mean given that distribution is CDF and sample type is continuous\n
         * <X> \approx \sum_bins (1.0 - P(X < bin_center)) * bin_width \n
         * yet to be implemented
         */
    }

    /** 
     * print mean of population probability distributions versus time \n
     * better to calculate moments directly using classes derived from Moment_Statistics\n
     */
    template<class Histogram_type, class Configuration_type>
    const double print_mean_populations(
    const Distribution_Statistics_Population<Histogram_type, Configuration_type> &statistics,
    const Pop &pop,
    const Spe &spe) {

        typedef typename Distribution_Statistics_Population<Histogram_type, Configuration_type>::Results_t Results_t;
        const array::Array3D<Results_t> array3D_results = statistics.get_results_so_far();

        std::string filename =
                "meanPopulationVersusTime__pop" + boost::lexical_cast<std::string > (pop.value()) +
                "__spe" + boost::lexical_cast<std::string > (spe.value()) +
                ".dat";
        boost::shared_ptr<std::ofstream> ofstream_ptr = monte_carlo::open_file_for_output(filename);

        for (int node = 0; node < array3D_results.get_dim2(); node++) {

            typedef typename Distribution_Statistics_Population<Histogram_type, Configuration_type>::Configuration_t::time_t time_type;
            const time_type tt = statistics.time_grid().at(Node(node));
            *ofstream_ptr << std::setw(30) << std::setprecision(20) << tt;

            typename Distribution_Statistics_Population<Histogram_type, Configuration_type>::Histogram_t::category distribution_category;
            const double mean = mean_detail::calculate_mean(
                    array3D_results.at(pop.value(), spe.value(), node),
                    distribution_category);
            *ofstream_ptr << std::setw(30) << std::setprecision(20) << mean;

            *ofstream_ptr << std::endl;
        }
    }
}

#endif	/* MEAN_POP_DIST_H */

