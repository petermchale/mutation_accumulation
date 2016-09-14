#ifndef MOMENT_STATISTICS_POPULATION_H
#define	MOMENT_STATISTICS_POPULATION_H

#include <mutation_accumulation/configuration/utilities/population.h> // monte_carlo::calculate_population
#include <mutation_accumulation/configuration/configuration/time_grid.h> // monte_carlo::Uniform_Time_Grid

#include "moment_statistics.h"

/*************************************************************************/

namespace monte_carlo {

    /** 
     * calculate indicated moment of population size of Spe-mutant in Pop sub-population at time indexed by Node 
     */
    template<class Moment_type, class Configuration_type>
    class Moment_Statistics_Population : public Moment_Statistics<Moment_type, Configuration_type> {
    private:

        typedef typename Configuration_type::population_t population_type;
        typedef typename Configuration_type::time_t time_type;

        typedef Moment_Statistics<Moment_type, Configuration_type> base_type;

    private:

        Uniform_Time_Grid<time_type> _time_grid;

    public:

        /**
         * constructor \n
         */
        explicit Moment_Statistics_Population(
                const Number_Pop &number_pop,
                const Number_Spe &number_spe,
                const Uniform_Time_Grid<time_type> &time_grid_)
        : base_type(
        number_pop,
        number_spe.value(),
        time_grid_.size(),
        number_spe.value(),
        time_grid_.size(),
        Pop(0),
        0,
        0,
        1e3), _time_grid(time_grid_) {

            /* sample type of Moment_type and population type of Configuration_type should agree */
            BOOST_STATIC_ASSERT((boost::is_same<typename Moment_type::sample_t, typename Configuration_type::population_t>::value));

        }

        /** 
         * dump results of a particular trial
         */
        virtual void dump(const Configuration_type &configuration) {

            for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
                for (int spe = 0; spe < configuration.number_species(); spe++)
                    for (int node = 0; node < configuration.number_nodes(); node++) {
                        const population_type population = calculate_population(configuration, Pop(pop), Spe(spe), Node(node));
                        this->update_moments(Pop(pop), spe, node, population);
                    }


            for (int spe = 0; spe < configuration.number_species(); spe++)
                for (int node = 0; node < configuration.number_nodes(); node++) {
                    const population_type population = calculate_population(configuration, Spe(spe), Node(node));
                    this->update_moments_whole(spe, node, population);
                }


        }

        const Uniform_Time_Grid<time_type> time_grid() const {

            return _time_grid;
        }

    };



}

#endif	/* MOMENT_STATISTICS_POPULATION_H */

