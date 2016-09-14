#ifndef FATE_H
#define	FATE_H

#include <mutation_accumulation/parameters/parameters_fwd.h> // Pop, etc
#include <mutation_accumulation/configuration/utilities/fpt_mutation.h> // mutation_occurred
#include <mutation_accumulation/configuration/utilities/fpt_extinction.h> // extinguished

/*************************************************************************/

namespace monte_carlo {

    /**
     * returns 1 if species spe arose in sub-population pop\n
     * returns 0 if species did not arise AND sub-population extinguished\n
     * return -1 if neither fate occurred
     * only tested for spe = 1 (starting from one type-0 stem cell per sub-population) 
     */
    template <class configuration_type>
    const int fate(
    const configuration_type &configuration,
    const Pop &pop,
    const Spe &spe) {

        if (mutation_occurred(configuration, pop, spe)) {
            return 1;
        } else if (!mutation_occurred(configuration, pop, spe) && extinguished(configuration, pop)) {
            return 0;
        } else {
            return -1;
        }
    }

    /**
     * returns 1 if species spe arose in any sub-population\n
     * returns 0 if species did not arise in any sub-population AND entire population extinguished\n
     * return -1 if neither fate occurred
     * only tested for spe = 1 (starting from one type-0 stem cell in each sub-population) 
     */
    template <class configuration_type>
    const int fate(
    const configuration_type &configuration,
    const Spe &spe) {

        if (mutation_occurred(configuration, spe)) {
            return 1;
        } else if (!mutation_occurred(configuration, spe) && extinguished(configuration)) {
            return 0;
        } else {
            return -1;
        }
    }

    /**
     * returns true if species spe arose in sub-population pop\n
     * returns false if species did not arise AND sub-population extinguished\n
     * produces error otherwise
     */
    template <class configuration_type>
    const bool fate_bool(
    const configuration_type &configuration,
    const Pop &pop,
    const Spe &spe) {

        if (fate(configuration, pop, spe) == 1) {
            return true;
        } else if (fate(configuration, pop, spe) == 0) {
            return false;
        } else if (fate(configuration, pop, spe) == -1) {
            std::cerr << "pop " << pop.value() << " has not acquired spe " << spe.value() << " nor has it extinguished" << std::endl;
            std::cerr << "configuration = " << std::endl;
            std::cerr << configuration << std::endl;
            assert(false);
        } else {
            std::cerr << "fate has returned invalid value" << std::endl;
            assert(false);
        }
    }

    /**
     * returns true if species spe arose in any sub-population\n
     * returns false if species did not arise in any sub-population AND entire population extinguished\n
     * produces error otherwise
     */
    template <class configuration_type>
    const bool fate_bool(
    const configuration_type &configuration,
    const Spe &spe) {

        if (fate(configuration, spe) == 1) {
            return true;
        } else if (fate(configuration, spe) == 0) {
            return false;
        } else if (fate(configuration, spe) == -1) {
            std::cerr << "entire population has not acquired spe " << spe.value() << " nor has it extinguished" << std::endl;
            std::cerr << "configuration = " << std::endl;
            std::cerr << configuration << std::endl;
            assert(false);
        } else {
            std::cerr << "fate has returned invalid value" << std::endl;
            assert(false);
        }
    }

    /**
     * returns true if one of two possible final fates (relative to last species) has occurred in sub-population pop\n
     * returns false otherwise\n
     */
    template <class configuration_type>
    const bool eventualFateOccurred_lastSpecies(
    const configuration_type &configuration,
    const Pop &pop) {

        const Spe last_species(configuration.number_species() - 1);

        const bool neither_fate_has_occurred = fate(configuration, pop, last_species) == -1;
        if (neither_fate_has_occurred) {
            return false;
        } else {
            return true;
        }
    }

    /**
     * returns true if final fate (relative to last species) has occurred in all sub-populations\n
     */
    template <class configuration_type>
    const bool eventualFateOccurred_lastSpecies_allSubPops(
    const configuration_type &configuration) {

        int final_fate_has_occurred_in_all_pops = 1;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
            final_fate_has_occurred_in_all_pops *= eventualFateOccurred_lastSpecies(configuration, Pop(pop));

        return final_fate_has_occurred_in_all_pops;
    }

}

#endif	/* FATE_H */

