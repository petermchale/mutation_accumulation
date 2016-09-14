#ifndef STATISTICS_EXTINCTION_H
#define	STATISTICS_EXTINCTION_H

#include <mutation_accumulation/probability/sample_space.h> // probability::make_uniform_sample_space
#include <mutation_accumulation/configuration/utilities/fpt_extinction.h>  // monte_carlo::extinction_time_whole

#include "distribution_statistics.h"

/*************************************************************************/

namespace monte_carlo {

    /** 
     * calculate distribution of time at which pop sub-population goes extinct
     */
    template<class Histogram_type, class Configuration_type>
    class Statistics_Extinction : public Distribution_Statistics<Histogram_type, Configuration_type> {
    private:

        typedef typename Configuration_type::time_t time_type; // Item 42
        typedef Distribution_Statistics<Histogram_type, Configuration_type> base_type;

    public:

        /**
         * constructor \n
         */
        explicit Statistics_Extinction(
                const Number_Pop &number_pop,
                const time_type &time_span,
                const double &error_probability,
                const Pop &pop_to_observe,
                const int &observer_divisor)
        : base_type(
        number_pop,
        1,
        1,
        1,
        1,
        probability::make_uniform_sample_space(51, time_span),
        error_probability,
        pop_to_observe,
        0,
        0,
        observer_divisor,
        "xxx",
        "xxx") {

            /* sample type of Histogram_type and time type of Configuration_type should agree */
            BOOST_STATIC_ASSERT((boost::is_same<typename Histogram_type::sample_t, typename Configuration_type::time_t>::value));

        }

        /** 
         * dump results of a particular trial
         */
        virtual void dump(const Configuration_type &configuration) {

            for (int pop = 0; pop < configuration.number_sub_pops(); pop++) {
                const time_type extinction_time = configuration.extinction_times(Pop(pop));
                this->update_histograms(Pop(pop), 0, 0, extinction_time);
            }

            {
                const time_type extinction_time = monte_carlo::extinction_time_whole(configuration);
                this->update_histograms_whole(0, 0, extinction_time);
            }

        }


    };



}

#endif	/* STATISTICS_EXTINCTION_H */

