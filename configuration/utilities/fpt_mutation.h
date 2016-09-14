#ifndef FPT_MUTATION_H
#define	FPT_MUTATION_H

#include <vector> // std::vector
#include <iomanip> // std::setw, etc

#include <mutation_accumulation/parameters/parameters_fwd.h> // Pop, etc

#include "lifetime_risk.h" // monte_carlo::mutation_occurred_within_timeSpan

/*************************************************************************/

namespace monte_carlo {

    /**
     * returns true if species spe has arisen in sub-population pop\n
     */
    template <class configuration_type>
    const bool mutation_occurred(
    const configuration_type &configuration,
    const Pop &pop,
    const Spe &spe) {

        typedef typename configuration_type::time_t time_type;

        /* convoluted implementation insures that a zero mutation time returns true \n
         * should implement this as T >= 0 (discrete) and T > -1e-6 (cts) AND test new implementation */
        return !(configuration.mutation_times(pop, spe) < static_cast<time_type> (0));
    }

    /**
     * returns true if species spe has arisen in at least one sub-population \n
     */
    template <class configuration_type>
    const bool mutation_occurred(
    const configuration_type &configuration,
    const Spe &spe) {

        int mutation_has_not_occurred = 1;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
            mutation_has_not_occurred *= !mutation_occurred(configuration, Pop(pop), spe);

        return !mutation_has_not_occurred;
    }

    /**
     * returns true if species spe has arisen in all sub-populations \n
     */
    template <class configuration_type>
    const bool all_mutations_occurred(
    const configuration_type &configuration,
    const Spe &spe) {

        int all_mutations_have_occurred = 1;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
            all_mutations_have_occurred *= mutation_occurred(configuration, Pop(pop), spe);

        return all_mutations_have_occurred;
    }

    /**
     * new: returns true if species spe has arisen in all sub-populations within time span defined by Time Grid member of Configuration\n
     */
    template <class configuration_type>
    const bool all_mutations_occurred__within_timeSpan(
    const configuration_type &configuration,
    const Spe &spe) {

        int all_mutations_have_occurred__within_timeSpan = 1;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
            all_mutations_have_occurred__within_timeSpan *= mutation_occurred_within_timeSpan_lessStringent(configuration, Pop(pop), spe);

        return all_mutations_have_occurred__within_timeSpan;
    }

    /**
     * time at which first mutant with spe mutations appeared in any sub-population
     */
    template <class configuration_type>
    const typename configuration_type::time_t mutation_time_whole(
    const configuration_type &configuration,
    const Spe &spe) {

        typedef typename configuration_type::time_t time_type;

        if (mutation_occurred(configuration, spe)) { // ... then mutation time is the smallest non-negative mutation time 

            // find non-negative mutation times
            std::vector<time_type> mutation_times; // empty vector
            for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
                if (mutation_occurred(configuration, Pop(pop), Spe(spe)))
                    mutation_times.push_back(configuration.mutation_times(Pop(pop), Spe(spe)));

            // find smallest non-negative mutation time
            const time_type very_large_time = static_cast<time_type>(1e8);
            time_type mutate_time_whole = very_large_time;
            for (int pop = 0; pop < mutation_times.size(); pop++) {
                const time_type mutation_time = mutation_times.at(pop);
                assert(!(mutation_time < static_cast<time_type> (0))); // all mutation times should be zero or greater
                if (mutation_time < mutate_time_whole)
                    mutate_time_whole = mutation_time;
            }

            return mutate_time_whole;

        } else { // ... then mutation time is -1

            return static_cast<time_type> (-1);
        }

    }

    /**
     * store mutation times to disk
     */
    template <class configuration_type>
    void store_mutation_times(
    configuration_type &configuration,
    const std::string &fileName) {

        std::ofstream fout(fileName.c_str()); // file closed when fout goes out of scope (RAII)

        if (!fout) {
            std::cerr << "cannot open " << fileName << std::endl;
            assert(fout);
        }

        for (int pop = 0; pop < configuration.number_sub_pops(); pop++) {

            for (int spe = 0; spe < configuration.number_species(); spe++)
                fout << std::setw(20) << std::setprecision(10) << configuration.mutation_times(Pop(pop), Spe(spe));

            fout << std::endl;
        }

    }


}

#endif	/* FPT_MUTATION_H */

