#ifndef CREATE_STATISTICS_H
#define	CREATE_STATISTICS_H

#include <mutation_accumulation/simulation/read_policy.h> 

/*************************************************************************/

namespace monte_carlo {

    /** 
     * create instance of a statistics gatherer 
     */
    template <class Statistics_Policy, class Configuration_Policy>
    const Statistics_Policy create_statistics(
            const Read_Policy_Base<Configuration_Policy> &read_policy) {

        return Statistics_Policy(
                Number_Pop(read_policy.get_population().number_sub_pops()),
                Number_Spe(read_policy.get_population().number_species()),
                read_policy.getTime_span_histogram(),
                read_policy.get_error_prob(),
                Pop(0), // observe statistics in first sub-population
                Spe(read_policy.get_population().number_species() - 1), // observe statistics in last species
                read_policy.get_observer_divisor());


    }

}

#endif	/* CREATE_STATISTICS_H */

