/* 
 * File:   cdf.h
 * Author: petermchale
 *
 * Created on June 5, 2012, 9:22 AM
 */

#ifndef CDF_H
#define	CDF_H

#include <mutation_accumulation/utility/data_traits.h> // data_types::discrete_type
#include <mutation_accumulation/utility/distribution_traits.h> // distribution_types::cdf_type

#include "histogram.h"

/*************************************************************************/

namespace probability {

    /**
     * cumulative probability distribution \n     
     * CDF(x) = P(X<=x)\n
     * code assumes X >= 0
     */
    template <class Notification_Policy>
    class CDF : public Histogram<typename Notification_Policy::sample_t> {
    public:

        typedef distribution_types::cdf_type category;

    private:

        typedef typename Notification_Policy::sample_t sample_type;
        typedef Histogram<sample_type> base_type;

    private:

        /**
         *  determine if event occurred given that sample is discrete
         */
        const bool event_occurred(const sample_type &sample, const sample_type &xx, data_types::discrete_type) {

            const bool cond1 = sample <= xx;
            const bool cond2 = sample >= 0;

            return (cond1 && cond2);
        }

        /**
         *  determine if event occurred given that sample is continuous
         */
        const bool event_occurred(const sample_type &sample, const sample_type &xx, data_types::continuous_type) {

            const bool cond1 = sample < xx;
            const bool cond2 = sample > -1e-8; // insures that a zero mutation time returns true

            return (cond1 && cond2);
        }

        /**
         * update CDF and notify observers
         */
        virtual void updateHistogram_notify(const sample_type &sample) {

            /* update CDF */
            for (int ii = 0; ii < this->get_sample_space().size(); ii++) {
                const sample_type xx = this->get_sample_space().at(ii);
                typename data_types::data_traits<sample_type>::category sample_category;
                if (event_occurred(sample, xx, sample_category))
                    this->increment_histogram(ii);
            }

            /* notify observers */
            if (Notification_Policy::notify(sample, this->get_sample_space()))
                this->notify();

        }

    public:

        /**
         * constructor
         */
        explicit CDF(const std::vector<sample_type> &sample_space_, const std::string &fileName) : base_type(sample_space_, fileName) {

        }

        /**
         * default constructor
         */
        explicit CDF() : base_type() {

        }
    };

}


#endif	

