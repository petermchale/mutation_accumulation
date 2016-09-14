#ifndef NOTIFICATION_POLICY_H
#define	NOTIFICATION_POLICY_H

#include <cassert> // assert 

#include <mutation_accumulation/utility/data_traits.h> // data_types::discrete_type, etc
#include <mutation_accumulation/configuration/utilities/lifetime_risk.h> // time_greater_than_zero

/*************************************************************************/

namespace probability {

    /**
     * rules dictating when probability subjects should notify observers\n
     * using independent classes rather than a hierarchy of classes because I want to invoke without instantiation
     */

    namespace notify_detail {

        /** 
         * implementation of Notify_Range::notify(..) for discrete sample type 
         */
        template <class sample_type>
        const bool do_range_notify(
        const sample_type &sample,
        const std::vector<sample_type> &sample_space,
        data_types::discrete_type) {

            const bool cond1 = sample >= sample_space.front();
            const bool cond2 = sample <= sample_space.back();

            return (cond1 && cond2);
        }

        /** 
         * implementation of Notify_Range::notify(..) for continuous sample type 
         */
        template <class sample_type>
        const bool do_range_notify(
        const sample_type &sample,
        const std::vector<sample_type> &sample_space,
        data_types::continuous_type) {

            const bool cond1 = sample > sample_space.front() - 1e-8; // insures that sample = sample_space.front() returns true, eg. sample = 0
            const bool cond2 = sample < sample_space.back();

            return (cond1 && cond2);
        }
    }

    /** 
     * notify observers if sample lies anywhere in given range \n
     */
    template <class sample_type>
    class Notify_Range {
    public:

        typedef sample_type sample_t;

        static const bool notify(const sample_type &sample, const std::vector<sample_type> &sample_space) {

            typename data_types::data_traits<sample_type>::category sample_category;

            notify_detail::do_range_notify(sample, sample_space, sample_category);

        }

    };

    /** 
     * notify observers if sample is true
     */
    class Notify_True {
    public:

        typedef bool sample_t;

        static const bool notify(const sample_t &sample, const std::vector<sample_t> &sample_space) {

            /* check that sample space is (0,1) */
            assert(sample_space.size() == 2);
            assert(sample_space.at(0) == false);
            assert(sample_space.at(1) == true);

            return (sample == true);

        }

    };

    /** 
     * notify observers if sample is non-negative and no larger than upper limit of sample space \n
     */
    template <class sample_type>
    class Notify_NonNegative_BoundedAbove {
    public:

        typedef sample_type sample_t;

        static const bool notify(const sample_type &sample, const std::vector<sample_type> &sample_space) {

            typename data_types::data_traits<sample_type>::category sample_category;
            const bool cond_lower = monte_carlo::mutation_occurred_within_timeSpan_detail::time_greater_than_zero(sample, sample_category);
            const bool cond_upper = monte_carlo::mutation_occurred_within_timeSpan_detail::time_less_than_LL(sample, sample_space.back(), sample_category);
            return cond_lower && cond_upper;

        }

    };


}

#endif	/* NOTIFICATION_POLICY_H */

