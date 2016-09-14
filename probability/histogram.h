#ifndef HISTOGRAM_BASE_H
#define	HISTOGRAM_BASE_H

#include <iomanip> // std::setw, etc
#include <algorithm> // std::max_element
#include <fstream> // std::ofstream
#include <iostream> // std::cerr
#include <cassert> // assert

#include <boost/date_time/posix_time/posix_time.hpp> // boost::posix_time

#include <mutation_accumulation/patterns/observer.h> // patterns::Subject
#include <mutation_accumulation/simulation/files.h> // monte_carlo::open_file_for_output 

/*************************************************************************/

/* consider moving implementation of large, implicitly inlined member functions to implementation file
 * (Item 30) */


/**
 * this namespace includes classes that empirically calculate various probability distributions
 */
namespace probability {

    /**
     * convenience class that contains a sample space and corresponding probability distribution\n
     * this class is copyable, unlike Histogram
     */
    template <class sample_type>
    class SampleSpace_Probability {
    public:

        std::vector<sample_type> _sample_space;
        std::vector<double> _probability;

        explicit SampleSpace_Probability(
                const std::vector<sample_type> &sample_space,
                const std::vector<double> &probability)
        :
        _sample_space(sample_space),
        _probability(probability) {

        }

        explicit SampleSpace_Probability() {

        }

        const std::vector<sample_type> sample_space() const {
            return _sample_space;
        }

        const std::vector<double> probability() const {
            return _probability;
        }

    };

    /**
     * abstract base class that holds the frequencies associated with a given sample space\n
     * cannot instantiate this class
     *
     * many of the member functions of this class could be made non-member (Item 23)\n
     */
    template <class sample_type>
    class Histogram : public patterns::Subject {
    public:

        typedef sample_type sample_t;
        typedef long long int number_trials_t;
        typedef long long int frequency_t;
        
    private:

        const std::vector<sample_type> sample_space;
        std::vector<frequency_t> histogram;
        std::vector<double> probability;
        number_trials_t number_trials;

        const std::string _fileName;

    private:

        typedef boost::posix_time::time_duration Duration_t;
        typedef boost::posix_time::ptime Posix_Time_t;
        typedef boost::posix_time::minutes Minutes_t;
        typedef boost::posix_time::second_clock Second_Clock_t;

        Duration_t interval; // interval between checkpoints
        Posix_Time_t next_checkpoint; // next checkpoint

    private:

        /* implementation member functions available to Histogram */

        /**
         * normalize histogram
         */
        void normalize() {

            if (number_trials > static_cast<number_trials_t>(0)) {

                for (int ii = 0; ii < sample_space.size(); ii++)
                    probability.at(ii) = (double) histogram.at(ii) / (double) number_trials;
            }

        }

        /**
         * write probability to disk
         */
        void store() {

            /* open file for (over-)writing */
            const boost::shared_ptr<std::ofstream> ofstream_ptr = monte_carlo::open_file_for_output(_fileName);

            for (int ii = 0; ii < probability.size(); ii++) {

                *ofstream_ptr << std::setw(10) << std::setprecision(3) << sample_space.at(ii);
                *ofstream_ptr << std::setw(20) << std::setprecision(10) << probability.at(ii);
                *ofstream_ptr << std::endl;

            }
        }

        /**
         * returns true if just passed a temporal checkpoint\n
         * should make a base class called Timer with virtual function onTick: Item 39\n
         * Histogram would contain a HistogramTimer class inheriting from Timer (Item 39)
         */
        const bool just_passed_checkpoint() {

            /* determine if time since start has exceeded current check point */
            const bool return_value = next_checkpoint < Second_Clock_t::local_time();

            if (return_value == true) {

                /* increment check point until it is greater than current time */
                while (next_checkpoint < Second_Clock_t::local_time())
                    next_checkpoint += interval;

            }

            return return_value;

        }


    protected:

        /* implementation member functions available to Histogram and derived classes */

        /**
         * constructor
         */
        explicit Histogram(const std::vector<sample_type> &sample_space_, const std::string &fileName)
        : patterns::Subject(), sample_space(sample_space_), number_trials(static_cast<number_trials_t>(0)), _fileName(fileName) {

            histogram = std::vector<frequency_t>(sample_space_.size(), 0);
            probability = std::vector<double>(sample_space_.size(), -1.0);

            interval = Minutes_t(15);
            next_checkpoint = Second_Clock_t::local_time() + interval;

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist\n
         *
         * this is an empty constructor\n
         * compiler should generate code to construct private data members (Item 30)
         */
        explicit Histogram() {

        }

        /**
         * increment histogram at position ii\n
         */
        void increment_histogram(const int &ii) {

            /* increment histogram */
            histogram.at(ii)++;

        }


        /**
         * update histogram of frequencies and notify any observers
         */
        virtual void updateHistogram_notify(const sample_type &sample) = 0; // implement in derived classes


    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Histogram() {

            store();
        }

        /**
         * update\n
         *
         */
        void update(const sample_type &sample) {

            /* update internal state */
            {
                /* update number of trials */
                number_trials++;

                /* update histogram of frequencies and notify any observers */
                updateHistogram_notify(sample);

                /* update probabilities */
                normalize();
            }

            /* write current state to disk at regular intervals */
            if (just_passed_checkpoint())
                store();

        }

        /**
         * get largest frequency
         */
        const frequency_t get_largest_frequency() const {

            return *(std::max_element(histogram.begin(), histogram.end()));

        }

        /**
         * get end frequency
         */
        const frequency_t get_end_frequency() const {

            return histogram.back();

        }

        /**
         * get sample space
         */
        const std::vector<sample_type> get_sample_space() const {

            return sample_space;

        }

        /**
         * get probability
         */
        const std::vector<double> get_probability() const {

            return probability;

        }

        /**
         * return number trials
         */
        const number_trials_t get_number_trials() const {

            return number_trials;

        }

    };






} // end namespace probability

#endif	/* HISTOGRAM_BASE_H */


