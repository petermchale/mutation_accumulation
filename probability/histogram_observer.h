#ifndef HISTOGRAM_OBSERVER_H
#define	HISTOGRAM_OBSERVER_H

#include <iomanip> // std::setw, etc
#include <fstream> // std::ofstream

#include <mutation_accumulation/patterns/observer.h> // patterns::Observer
#include <mutation_accumulation/utility/distribution_traits.h> // distribution_types::

/*************************************************************************/

namespace probability {

    /**
     * Histogram_Observer is notified when a Histogram of type Histogram_type is changed \n
     */
    template <class Histogram_type>
    class Histogram_Observer : public patterns::Observer {
    private:

        /**
         * To set _histogram to "null" in eg. a default ctor, one needs a pointer, not a reference\n
         * A reference must be initialised to refer to something\n
         */
        Histogram_type &_histogram;

        std::ofstream flog; // log file

        int number_updates; // keep track of number of times update() has been called
        const int divisor; // print to log file every divisor updates

    private:

        /** 
         * print log data using end frequency \n
         * this member function is non-const because it modifies member variable flog
         */
        void print_log_data__end_frequency() {

            /* query subject */
            const typename Histogram_type::frequency_t end_frequency = _histogram.get_end_frequency();
            const typename Histogram_type::number_trials_t number_trials = _histogram.get_number_trials();

            flog << std::setw(20) << end_frequency;
            flog << std::setw(20) << number_trials;
            flog << std::setw(20) << (double) end_frequency / (double) number_trials;
            flog << std::endl;
        }

        /** 
         * print log data given that distribution is CDF \n
         * this member function is non-const because it modifies member variable flog
         */
        void print_log_data(distribution_types::cdf_type) {

            print_log_data__end_frequency();

        }

        /** 
         * print log data given that distribution is PMF with a boolean sample space\n
         * this member function is non-const because it modifies member variable flog\n
         */
        void print_log_data(distribution_types::pmf_bool_type) {

            print_log_data__end_frequency();

        }

    public:

        /**
         * constructor
         */
        explicit Histogram_Observer(Histogram_type &histogram, const std::string &logFileName, const int &_divisor)
        : patterns::Observer(), _histogram(histogram), number_updates(0), divisor(_divisor) {

            /* open a log file */
            flog.open(logFileName.c_str(), std::ios_base::out); // file closed when flog goes out of scope (RAII)

            if (!flog.is_open()) {
                std::cerr << "cannot open " << logFileName << std::endl;
                assert(false);
            }


        }

        /**
         * virtual destructor\n
         */
        virtual ~Histogram_Observer() {

            _histogram.detach(this);

        }

        /**
         * attach current object to the list of observers in subject \n
         * cannot call this in ctor because object must be constructed before "this" is called\n
         */
        void attach_to_subject() {

            _histogram.attach(this);

        }

        /**
         * update by querying subject\n
         */
        virtual void update() {

            /* update number of updates */
            number_updates++;

            /* write to log file every divisor updates */
            if (number_updates % divisor == 0) {

                typename Histogram_type::category distribution_category;
                print_log_data(distribution_category);

            }

        }

    };



}


#endif	

