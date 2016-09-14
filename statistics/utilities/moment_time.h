#ifndef MOMENT_TIME_H
#define	MOMENT_TIME_H

#include <cassert> // assert 

/*************************************************************************/

namespace monte_carlo {

    namespace moment_time_detail { 
        
        template<class Statistics_Policy>
        void check_dim2(const Statistics_Policy &statistics) {
            
            assert(statistics.get_results_so_far().get_dim2() == 1);
        }
    }
    
    /** 
     * moment of mutation time for species spe in sub-population pop \n
     * implement similar function for Moment_Statistics_ExtinctionTime omitting Spe argument \n
     * check that dim1 == 1 \n
     * factor out code common to both functions
     */
    template<class Moment_type, class Configuration_type>
    const typename Moment_type::Results_t
    moment_time(
    const Moment_Statistics_MutationTime<Moment_type, Configuration_type> &statistics,
    const Pop &pop,
    const Spe &spe) {

        moment_time_detail::check_dim2(statistics);
        
        return statistics.get_results_so_far().at(pop.value(), spe.value(), 0);

    }

    /** 
     * moment of mutation time for species spe in whole population \n
     */
    template<class Moment_type, class Configuration_type>
    const typename Moment_type::Results_t
    moment_time(
    const Moment_Statistics_MutationTime<Moment_type, Configuration_type> &statistics,
    const Spe &spe) {

        moment_time_detail::check_dim2(statistics);

        return statistics.get_results_so_far_whole().at(spe.value(), 0);

    }
}

#endif	/* MOMENT_TIME_H */

