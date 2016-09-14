#ifndef POPULATION_H
#define	POPULATION_H

#include <cassert> // assert

/*************************************************************************/

namespace monte_carlo {

    /**
     * population size specific to pop, spe, and node for a complete path
     */
    template <class configuration_type>
    const typename configuration_type::population_t calculate_population(
    const configuration_type &configuration,
    const Pop &pop,
    const Spe &spe,
    const Node &node) {

        return configuration.path().at(pop, spe, node);
    }

    /**
     * population size specific to spe and node for a complete path
     */
    template <class configuration_type>
    const typename configuration_type::population_t calculate_population(
    const configuration_type &configuration,
    const Spe &spe,
    const Node &node) {

        typedef typename configuration_type::population_t population_type;

        population_type population = static_cast<population_type> (0);
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
            population += configuration.path().at(Pop(pop), spe, node);

        return population;
    }
    
    
}

#endif	/* POPULATION_H */

