#ifndef STATISTICS_MUTATION_H
#define	STATISTICS_MUTATION_H

#include <mutation_accumulation/probability/sample_space.h> // probability::make_uniform_sample_space
#include <mutation_accumulation/configuration/utilities/fpt_mutation.h>  // monte_carlo::mutation_time_whole

#include "distribution_statistics.h" 

/*************************************************************************/

namespace monte_carlo {

    /** 
     * calculate distributions of first time at which spe-mutant arises in pop sub-population 
     */
    template<class Histogram_type, class Configuration_type>
    class Statistics_Mutation : public Distribution_Statistics<Histogram_type, Configuration_type> {
    private:

        typedef typename Configuration_type::time_t time_type; // Item 42

        typedef Distribution_Statistics<Histogram_type, Configuration_type> base_type;

    private:

        void error_checking() const {

            /* sample type of Histogram_type and time type of Configuration_type should agree */
            BOOST_STATIC_ASSERT((boost::is_same<typename Histogram_type::sample_t, typename Configuration_type::time_t>::value));

        }
    public:

        /**
         * constructor \n
         */
        explicit Statistics_Mutation(
                const Number_Pop &number_pop,
                const Number_Spe &number_spe,
                const time_type &sample_space_span,
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
        grid::make_uniform_grid(100, static_cast<time_type> (0), sample_space_span),
        error_probability,
        pop_to_observe,
        spe_to_observe.value(),
        0,
        observer_divisor,
        "spe",
        "xxx") {

            error_checking();

        }

        /**
         * constructor \n
         */
        explicit Statistics_Mutation(
                const Number_Pop &number_pop,
                const Number_Spe &number_spe,
                const time_type &sample_space_lower,
                const time_type &sample_space_upper,
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
        grid::make_logarithmic_grid(100, sample_space_lower, sample_space_upper),
        error_probability,
        pop_to_observe,
        spe_to_observe.value(),
        0,
        observer_divisor,
        "spe",
        "xxx") {

            error_checking();

        }

        /** 
         * dump results of a particular trial
         */
        virtual void dump(const Configuration_type &configuration) {

            for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
                for (int spe = 0; spe < configuration.number_species(); spe++) {
                    const time_type mutation_time = configuration.mutation_times(Pop(pop), Spe(spe));                    
                    this->update_histograms(Pop(pop), spe, 0, mutation_time);
                }

            for (int spe = 0; spe < configuration.number_species(); spe++) {
                const time_type mutation_time = monte_carlo::mutation_time_whole(configuration, Spe(spe));
                this->update_histograms_whole(spe, 0, mutation_time);
            }


        }

    };

}

#endif	/* STATISTICS_MUTATION_H */

