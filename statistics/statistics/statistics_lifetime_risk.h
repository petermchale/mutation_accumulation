#ifndef STATISTICS_LIFETIME_RISK_H
#define	STATISTICS_LIFETIME_RISK_H

#include <mutation_accumulation/probability/sample_space.h> // probability::make_bernoulli_sample_space
#include <mutation_accumulation/configuration/utilities/lifetime_risk.h> // monte_carlo::mutation_occurred_within_timeSpan 
#include <mutation_accumulation/probability/pmf.h> 

#include "distribution_statistics.h" 

/*************************************************************************/

namespace monte_carlo {

    namespace Statistics_Lifetime_Risk_namespace {

        typedef probability::PMF_Bool Histogram_type;

        /** 
         * calculate the complementary probabilities that sub-population pop \n
         * accumulates spe mutations (sample = true) or doesn't (sample = false) \n
         * within time span defined by Time Grid member of Configuration \n
         * 
         */
        template<class Configuration_type>
        class Statistics_Lifetime_Risk : public Distribution_Statistics<Histogram_type, Configuration_type> {
        private:

            typedef Distribution_Statistics<Histogram_type, Configuration_type> base_type;

        public:

            /**
             * constructor \n
             */
            explicit Statistics_Lifetime_Risk(
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
                        const bool fate = monte_carlo::mutation_occurred_within_timeSpan(configuration, Pop(pop), Spe(spe));
                        this->update_histograms(Pop(pop), spe, 0, fate);
                    }

                for (int spe = 0; spe < configuration.number_species(); spe++) {
                    const bool fate = monte_carlo::mutation_occurred_within_timeSpan(configuration, Spe(spe));
                    this->update_histograms_whole(spe, 0, fate);
                }


            }

        };

    }

    using Statistics_Lifetime_Risk_namespace::Statistics_Lifetime_Risk;

    /** 
     * probability that sub-population pop accumulates spe mutations within time span (sample = true)
     */
    template<class Configuration_type>
    const double probability_mutation_fate(
    const Statistics_Lifetime_Risk<Configuration_type> &statistics,
    const Pop &pop,
    const Spe &spe) {

        return statistics.get_results_so_far().at(pop.value(), spe.value(), 0).probability().back();

    }

    /** 
     * probability that whole population accumulates spe mutations within time span (sample = true) \n
     */
    template<class Configuration_type>
    const double probability_mutation_fate(
    const Statistics_Lifetime_Risk<Configuration_type> &statistics,
    const Spe &spe) {

        return statistics.get_results_so_far_whole().at(spe.value(), 0).probability().back();

    }

    /** 
     * probability that whole population accumulates full complement of mutations within time span (sample = true) \n
     */
    template<class Configuration_type>
    const double probability_mutation_fate(
    const Statistics_Lifetime_Risk<Configuration_type> &statistics) {

        const Spe last_species(statistics.get_results_so_far_whole().get_dim0() - 1);

        return probability_mutation_fate(statistics, last_species);

    }

}

#endif	/* STATISTICS_LIFETIME_RISK_H */

