#ifndef FPT_EXTINCTION_H
#define	FPT_EXTINCTION_H

#include <mutation_accumulation/parameters/parameters_fwd.h> // Pop, etc

/*************************************************************************/

namespace monte_carlo {

    /**
     * returns true if sub-population pop has extinguished
     */
    template <class configuration_type>
    const bool extinguished(
    const configuration_type &configuration,
    const Pop &pop) {

        typedef typename configuration_type::time_t time_type;

        return (configuration.extinction_times(Pop(pop)) > static_cast<time_type> (0));

    }

    /**
     * returns true if all sub-populations have extinguished
     */
    template <class configuration_type>
    const bool extinguished(
    const configuration_type &configuration) {

        int all_pops_extinguished = 1;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++) {
            all_pops_extinguished *= extinguished(configuration, Pop(pop));
        }

        return all_pops_extinguished;
    }

    /**
     * time at which whole population extinguished 
     */
    template <class configuration_type>
    const typename configuration_type::time_t extinction_time_whole(const configuration_type &configuration) {
        
        typedef typename configuration_type::time_t time_type;
        
        if (extinguished(configuration)) {
            
            // extinction time is time at which last surviving sub-population extinguished
            time_type extinct_time_whole = -1;
            for (int pop = 0; pop < configuration.number_sub_pops(); pop++) {
                const time_type extinction_time = configuration.extinction_times(Pop(pop));
                if (extinction_time > extinct_time_whole) extinct_time_whole = extinction_time;
            }
            
            return extinct_time_whole;
            
        } else {
            
            return static_cast<time_type> (-1);
        }


    }


}

#endif	/* FPT_EXTINCTION_H */

