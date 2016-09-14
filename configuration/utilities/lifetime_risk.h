#ifndef CONVENIENCE_LIFETIME_RISK_H
#define	CONVENIENCE_LIFETIME_RISK_H

#include <mutation_accumulation/parameters/parameters_fwd.h> // Pop, etc
#include <mutation_accumulation/configuration/utilities/fpt_mutation.h> // monte_carlo::mutation_time_whole
#include <mutation_accumulation/configuration/utilities/print.h> // <<configuration

/*************************************************************************/

namespace monte_carlo {

    namespace mutation_occurred_within_timeSpan_detail {

        /** 
         * compare time to zero \n
         * discrete-time implementation \n
         * \n
         * this uses traits instead of simply matching f(int) versus f(double)\n
         * otherwise f(double) will call f(int) by type conversion if user forgets to define f(double)
         */
        template <class time_type>
        const bool time_greater_than_zero(const time_type &tt, data_types::discrete_type) {

            return (tt >= 0);
        }

        /** 
         * compare time to zero \n
         * continuous-time implementation \n
         * \n
         * this uses traits instead of simply matching f(int) versus f(double)\n
         * otherwise f(double) will call f(int) by type conversion if user forgets to define f(double)
         */
        template <class time_type>
        const bool time_greater_than_zero(const time_type &tt, data_types::continuous_type) {

            return (tt >= -1e-8);
        }

        /** 
         * compare time to life span \n
         * discrete-time implementation \n
         * \n
         * this uses traits instead of simply matching f(int) versus f(double)\n
         * otherwise f(double) will call f(int) by type conversion if user forgets to define f(double)
         */
        template <class time_type>
        const bool time_less_than_LL(const time_type &tt, const time_type &LL, data_types::discrete_type) {

            return (tt <= LL);
        }

        /** 
         * compare time to life span \n
         * cts-time implementation \n
         * \n
         * this uses traits instead of simply matching f(int) versus f(double)\n
         * otherwise f(double) will call f(int) by type conversion if user forgets to define f(double)
         */
        template <class time_type>
        const bool time_less_than_LL(const time_type &tt, const time_type &LL, data_types::continuous_type) {

            return (tt < LL);
        }

        /**
         * returns true if mutation_time lies in [0, time span defined by Time Grid member of Configuration] \n
         * this is intended to be called by a Statistics_Gatherer object AFTER a trajectory has been generated
         */
        template <class configuration_type>
        const bool do_mutation_occurred_within_timeSpan(
        const configuration_type &configuration,
        const typename configuration_type::time_t &mutation_time) {

            typename data_types::data_traits<typename configuration_type::time_t>::category time_category;
            if (time_greater_than_zero(mutation_time, time_category)) {

                if (time_less_than_LL(mutation_time, configuration.get_last_node_time(), time_category))
                    return true;
                else
                    return false;

            } else {

                if (configuration.end() || extinguished(configuration)) { // 8th March 2013
                    return false;
                } else {
                    std::cerr << "configuration = " << configuration << std::endl;
                    std::cerr << "mutation_time = " << mutation_time << std::endl;
                    std::cerr << "cannot determine if mutation occurred within time span because path was not complete" << std::endl;
                    assert(false);
                }

            }

        }

        /**
         * returns true if mutation_time lies in [0, time span defined by Time Grid member of Configuration] \n
         * this is intended to be called by a Path_Policy_Base object BEFORE AND AFTER a trajectory has been completed
         */
        template <class configuration_type>
        const bool do_mutation_occurred_within_timeSpan_lessStringent(
        const configuration_type &configuration,
        const typename configuration_type::time_t &mutation_time) {

            typename data_types::data_traits<typename configuration_type::time_t>::category time_category;
            if (time_greater_than_zero(mutation_time, time_category)) {

                if (time_less_than_LL(mutation_time, configuration.get_last_node_time(), time_category))
                    return true;
                else
                    return false;

            } else {

                return false;
            }

        }
    }

    /**
     * returns true if species spe has arisen in sub-population pop within time span defined by Time Grid member of Configuration \n
     */
    template <class configuration_type>
    const bool mutation_occurred_within_timeSpan(
    const configuration_type &configuration,
    const Pop &pop,
    const Spe &spe) {

        return mutation_occurred_within_timeSpan_detail::do_mutation_occurred_within_timeSpan(configuration, configuration.mutation_times(pop, spe));

    }

    /**
     * new: returns true if species spe has arisen in sub-population pop within time span defined by Time Grid member of Configuration \n
     */
    template <class configuration_type>
    const bool mutation_occurred_within_timeSpan_lessStringent(
    const configuration_type &configuration,
    const Pop &pop,
    const Spe &spe) {

        return mutation_occurred_within_timeSpan_detail::do_mutation_occurred_within_timeSpan_lessStringent(configuration, configuration.mutation_times(pop, spe));

    }

    /**
     * returns true if species spe has arisen in any sub-population within time span defined by Time Grid member of Configuration \n
     */
    template <class configuration_type>
    const bool mutation_occurred_within_timeSpan(
    const configuration_type &configuration,
    const Spe &spe) {

        return mutation_occurred_within_timeSpan_detail::do_mutation_occurred_within_timeSpan(configuration, mutation_time_whole(configuration, spe));

    }


}

#endif	/* CONVENIENCE_LIFETIME_RISK_H */

