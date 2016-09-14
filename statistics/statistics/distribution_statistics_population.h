#ifndef DISTRIBUTION_STATISTICS_POPULATION_H
#define	DISTRIBUTION_STATISTICS_POPULATION_H

#include <mutation_accumulation/probability/sample_space.h> // probability::make_uniform_sample_space
#include <mutation_accumulation/configuration/utilities/population.h> // monte_carlo::calculate_population
#include <mutation_accumulation/configuration/configuration/time_grid.h> // monte_carlo::Uniform_Time_Grid

#include "distribution_statistics.h"

/*************************************************************************/

namespace monte_carlo {

    /** 
     * calculate distribution of population size of Spe-mutant in Pop sub-population at time indexed by Node 
     */
    template<class Histogram_type, class Configuration_type>
    class Distribution_Statistics_Population : public Distribution_Statistics<Histogram_type, Configuration_type> {
    private:

        typedef typename Configuration_type::population_t population_type;
        typedef typename Configuration_type::time_t time_type;

        typedef Distribution_Statistics<Histogram_type, Configuration_type> base_type;

    private:

        Uniform_Time_Grid<time_type> _time_grid;

    public:

        /**
         * constructor \n
         */
        explicit Distribution_Statistics_Population(
                const Number_Pop &number_pop,
                const Number_Spe &number_spe,
                const Uniform_Time_Grid<time_type> &time_grid_,
                const population_type &population_span_histogram,
                const double &error_probability,
                const Pop &pop_to_observe,
                const Spe &spe_to_observe,
                const Node &node_to_observe,
                const int &observer_divisor)
        : base_type(
        number_pop,
        number_spe.value(),
        time_grid_.size(),
        number_spe.value(),
        time_grid_.size(),
        /* logarithmic space is necessary to plot interesting parts of histograms for each species simultaneously \n
         * sample space that includes all integral sample values is necessary to estimate mean from discrete CDF or PMF \n
         * "full" sample space is achieved by making suggested number of sample points larger than sample range */
        probability::make_uniform_sample_space(1e4, population_span_histogram),
        error_probability,
        pop_to_observe,
        spe_to_observe.value(),
        node_to_observe.value(),
        observer_divisor,
        "spe",
        "node"), _time_grid(time_grid_) {

            /* sample type of Histogram_type and population type of Configuration_type should agree */
            BOOST_STATIC_ASSERT((boost::is_same<typename Histogram_type::sample_t, typename Configuration_type::population_t>::value));

        }

        /** 
         * dump results of a particular trial
         */
        virtual void dump(const Configuration_type &configuration) {

            for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
                for (int spe = 0; spe < configuration.number_species(); spe++)
                    for (int node = 0; node < configuration.number_nodes(); node++) {
                        const population_type population = calculate_population(configuration, Pop(pop), Spe(spe), Node(node));
                        this->update_histograms(Pop(pop), spe, node, population);
                    }


            for (int spe = 0; spe < configuration.number_species(); spe++)
                for (int node = 0; node < configuration.number_nodes(); node++) {
                    const population_type population = calculate_population(configuration, Spe(spe), Node(node));
                    this->update_histograms_whole(spe, node, population);
                }


        }

        const Uniform_Time_Grid<time_type> time_grid() const {

            return _time_grid;
        }

    };



}

#endif	/* DISTRIBUTION_STATISTICS_POPULATION_H */

