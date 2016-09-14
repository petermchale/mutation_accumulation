#ifndef FPT_MUTATION_LASTSPECIES_H
#define	FPT_MUTATION_LASTSPECIES_H

#include <mutation_accumulation/parameters/parameters_fwd.h> // Pop, etc

#include "fpt_mutation.h" 

/*************************************************************************/

namespace monte_carlo {


    /**
     * time at which first mutant with full complement of mutations appeared in any sub-population
     */
    template <class configuration_type>
    const typename configuration_type::time_t lastSpecies_mutation_time(
    const configuration_type &configuration) {

        const Spe last_species(configuration.number_species() - 1);

        return mutation_time_whole(configuration, last_species);
    }

    /**
     * time at which first mutant with full complement of mutations appeared in sub-population pop
     */
    template <class configuration_type>
    const typename configuration_type::time_t lastSpecies_mutation_time(
    const configuration_type &configuration,
    const Pop &pop) {

        const Spe last_species(configuration.number_species() - 1);

        return configuration.mutation_times(pop, last_species);
    }

    /**
     * returns true if last species has arisen in at least one sub-population \n
     */
    template <class configuration_type>
    const bool last_species_occurred(
    const configuration_type &configuration) {

        const Spe last_species(configuration.number_species() - 1);

        return mutation_occurred(configuration, last_species);

    }

    /**
     * returns true if last species has arisen in all sub-populations \n
     */
    template <class configuration_type>
    const bool all_last_species_occurred(
    const configuration_type &configuration) {

        const Spe last_species(configuration.number_species() - 1);

        return all_mutations_occurred(configuration, last_species);

    }

    /**
     * new: returns true if last species has arisen in all sub-populations within time span defined by Time Grid member of Configuration \n
     */
    template <class configuration_type>
    const bool all_last_species_occurred__within_timeSpan(
    const configuration_type &configuration) {

        const Spe last_species(configuration.number_species() - 1);

        return all_mutations_occurred__within_timeSpan(configuration, last_species);

    }



}

#endif	/* FPT_MUTATION_LASTSPECIES_H */

