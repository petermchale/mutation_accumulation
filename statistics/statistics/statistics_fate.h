#ifndef STATISTICS_FATE_H
#define	STATISTICS_FATE_H

#include <mutation_accumulation/probability/sample_space.h> // probability::make_bernoulli_sample_space
#include <mutation_accumulation/configuration/utilities/fate.h> // monte_carlo::fate_bool 
#include <mutation_accumulation/probability/pmf.h> 

#include "distribution_statistics.h" 

/*************************************************************************/

namespace monte_carlo {

    namespace Statistics_Fate_namespace {

        typedef probability::PMF_Bool Histogram_type;

        /** 
         * calculate the complementary probabilities that sub-population pop \n
         * eventually accumulates spe mutations (sample = true) or \n
         * eventually extinguishes (sample = false)
         */
        template<class Configuration_type>
        class Statistics_Fate : public Distribution_Statistics<Histogram_type, Configuration_type> {
        private:

            typedef Distribution_Statistics<Histogram_type, Configuration_type> base_type;

        public:

            /**
             * constructor \n
             */
            explicit Statistics_Fate(
                    const Number_Pop &number_pop,
                    const Number_Spe &number_spe,
                    const double &error_probability,
                    const Pop &pop_to_observe,
                    const Spe &spe_to_observe,
                    const int &observer_divisor)
            : base_type(
            number_pop,
            number_spe.value(),
            1,
            number_spe.value(),
            1,
            probability::make_bernoulli_sample_space(),
            error_probability,
            pop_to_observe,
            spe_to_observe.value(),
            0,
            observer_divisor,
            "spe",
            "xxx") {

                /* sample type of Histogram_type should be bool */
                BOOST_STATIC_ASSERT((boost::is_same<typename Histogram_type::sample_t, bool>::value));

            }

            /** 
             * dump results of a particular trial
             */
            virtual void dump(const Configuration_type &configuration) {

                for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
                    for (int spe = 0; spe < configuration.number_species(); spe++) {
                        const bool fate = monte_carlo::fate_bool(configuration, Pop(pop), Spe(spe));
                        this->update_histograms(Pop(pop), spe, 0, fate);
                    }

                for (int spe = 0; spe < configuration.number_species(); spe++) {
                    const bool fate = monte_carlo::fate_bool(configuration, Spe(spe));
                    this->update_histograms_whole(spe, 0, fate);
                }


            }

        };

    }

    using Statistics_Fate_namespace::Statistics_Fate;

    /** 
     * probability that sub-population pop eventually accumulates spe mutations (sample = true)
     */
    template<class Configuration_type>
    const double probability_mutation_fate(
    const Statistics_Fate<Configuration_type> &statistics,
    const Pop &pop,
    const Spe &spe) {

        return statistics.get_results_so_far().at(pop.value(), spe.value(), 0).probability().back();

    }

    /** 
     * probability that whole population eventually accumulates spe mutations (sample = true) \n
     * statistics object contains a copy of a reference to original histogram objects
     */
    template<class Configuration_type>
    const double probability_mutation_fate(
    const Statistics_Fate<Configuration_type> &statistics,
    const Spe &spe) {

        return statistics.get_results_so_far_whole().at(spe.value(), 0).probability().back();

    }

    /** 
     * probability that whole population eventually gets full complement of mutations (sample = true) \n
     */
    template<class Configuration_type>
    const double probability_mutation_fate(
    const Statistics_Fate<Configuration_type> &statistics) {

        const Spe last_species(statistics.get_results_so_far_whole().get_dim0() - 1);

        return statistics.get_results_so_far_whole().at(last_species.value(), 0).probability().back();

    }

}

#endif	/* STATISTICS_FATE_H */

