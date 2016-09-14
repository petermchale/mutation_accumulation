#ifndef PMF_H
#define	PMF_H

//#define PMF_DEBUG

#include <mutation_accumulation/utility/data_traits.h> 
#include <mutation_accumulation/utility/distribution_traits.h> 

#include "histogram.h"
#include "notification_policy.h"

/*************************************************************************/

namespace probability {

    /**
     * class name: "PMF"
     * probability mass function: discrete sample space \n
     * PMF(x) = P(X is in bin centered on x)\n
     * not yet implemented (and will take some effort to implement it)
     * member function event_occurred (see below) should be altered to determine if sample lies in given bin
     */


    namespace PMF_Bool_namespace {

        typedef Notify_True Notification_Policy;

        /**
         * probability mass function: boolean sample space x = (false, true) \n
         * PMF(x) = P(X = x)\n
         */
        class PMF_Bool : public Histogram<Notification_Policy::sample_t> {
        public:

            typedef distribution_types::pmf_bool_type category;

        private:

            typedef Notification_Policy::sample_t sample_type;
            typedef Histogram<sample_type> base_type;

        private:

            /**
             * determine if event occurred on discrete sample space \n
             */
            const bool event_occurred(const sample_type &sample, const sample_type &xx, data_types::discrete_type) {

                return (sample == xx);
            }

            /**
             * update PMF and notify any observers
             */
            virtual void updateHistogram_notify(const sample_type &sample) {

#ifdef PMF_DEBUG 
                std::cout << "sample = " << sample << std::endl;
#endif            
                /* update PMF */
                for (int ii = 0; ii < this->get_sample_space().size(); ii++) {
                    const sample_type xx = this->get_sample_space().at(ii);
                    data_types::data_traits<sample_type>::category sample_category;
                    if (event_occurred(sample, xx, sample_category)) {
#ifdef PMF_DEBUG 
                        std::cout << "event occurred at " << xx << std::endl;
#endif
                        this->increment_histogram(ii);
                        break;
                    }
                }

                /* notify observers */
                if (Notification_Policy::notify(sample, this->get_sample_space())) {
#ifdef PMF_DEBUG 
                    std::cout << "notify: sample = " << sample << std::endl;
#endif 
                    this->notify();
                }

            }

        public:

            /**
             * constructor
             */
            explicit PMF_Bool(const std::vector<sample_type> &sample_space_, const std::string &fileName) : base_type(sample_space_, fileName) {

            }

            /**
             * default constructor
             */
            explicit PMF_Bool() : base_type() {

            }
        };
    }

    using PMF_Bool_namespace::PMF_Bool;
}


#endif	

