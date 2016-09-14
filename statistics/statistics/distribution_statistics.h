#ifndef DISTRIBUTION_STATISTICS_H
#define	DISTRIBUTION_STATISTICS_H

#include <boost/shared_ptr.hpp> // boost::shared_ptr 

#include <mutation_accumulation/probability/histogram_observer.h> // probability::Histogram_Observer
#include <mutation_accumulation/parameters/parameters_fwd.h> // monte_carlo::Number_Pop
#include <mutation_accumulation/patterns/observer.h> // patterns::createAttachedObserver
#include <mutation_accumulation/probability/histogram.h> // probability::SampleSpace_Probability

#include "statistics_gatherer.h" // monte_carlo::Statistics_Gatherer

/*************************************************************************/

namespace monte_carlo {

    /**
     * each sub-population is endowed with a matrix of histograms\n
     * should generalize matrix to array of arbitrary dimensionality\n
     */
    template<class Histogram_type, class Configuration_type>
    class Distribution_Statistics : public Statistics_Gatherer<Configuration_type, probability::SampleSpace_Probability<typename Histogram_type::sample_t> > {
    public:

        typedef Histogram_type Histogram_t;

    private:

        typedef boost::shared_ptr<Histogram_type> Histogram_ptr_type;
        typedef array::Array3D<Histogram_ptr_type> array3D_histograms_type;
        typedef array::Array2D<Histogram_ptr_type> array2D_histograms_type;

        typedef probability::Histogram_Observer<Histogram_type> Histogram_Observer_type;

        typedef typename Histogram_type::sample_t sample_type;
        typedef probability::SampleSpace_Probability<sample_type> SampleSpace_Probability_type;


    private:

        /**
         * pointers to histograms are necessary since histograms are uncopyable \n
         * histogram pointers are also used for consistency with Histogram_observer pointers returned by createAttachedObserver
         */
        array3D_histograms_type _histograms; // arbitrary matrix of histograms for each sub-population 
        array2D_histograms_type _histograms_whole; // independent arbitrary matrix of whole-population histograms

        boost::shared_ptr<Histogram_Observer_type> _histogram_observer; // observe a single histogram

        const double _error_probability; // error in probability 
        const int _pop_to_observe; // sub-population to monitor
        const int _jj_per_pop_to_observe; // which row of histogram matrix to monitor in sub-population _pop_to_observe 
        const int _kk_per_pop_to_observe; // which col of histogram matrix to monitor in sub-population _pop_to_observe 

    private:

        /** 
         * fetch sample space and probability distribution from histogram ptr
         */
        const SampleSpace_Probability_type fetch_sampleSpace_probability(const Histogram_ptr_type &histogram_ptr) const {

            const std::vector<sample_type> ss = histogram_ptr->get_sample_space();
            const std::vector<double> pp = histogram_ptr->get_probability();

            return SampleSpace_Probability_type(ss, pp);

        }


    protected:

        /**
         * constructor \n
         * base-class default ctor called by compiler
         */
        explicit Distribution_Statistics(
                const Number_Pop &number_pop,
                const int &dim1_per_pop,
                const int &dim2_per_pop,
                const int &dim0_whole,
                const int &dim1_whole,
                const std::vector<sample_type> &sample_space,
                const double &error_probability,
                const Pop &pop_to_observe,
                const int &jj_per_pop_to_observe,
                const int &kk_per_pop_to_observe,
                const int &observer_divisor,
                const std::string &label_for_dim1_per_pop,
                const std::string &label_for_dim2_per_pop)
        :
        _histograms(array3D_histograms_type(number_pop.value(), dim1_per_pop, dim2_per_pop)),
        _histograms_whole(array2D_histograms_type(dim0_whole, dim1_whole)),
        _error_probability(error_probability),
        _pop_to_observe(pop_to_observe.value()),
        _jj_per_pop_to_observe(jj_per_pop_to_observe),
        _kk_per_pop_to_observe(kk_per_pop_to_observe) {

            /* check template parameter types to supplement "duck typing" */
            BOOST_STATIC_ASSERT((boost::is_base_of<probability::Histogram<sample_type>, Histogram_type>::value));

            /* create sub-population histogram matrices */
            for (int pop = 0; pop < number_pop.value(); pop++)
                for (int jj = 0; jj < dim1_per_pop; jj++)
                    for (int kk = 0; kk < dim2_per_pop; kk++) {

                        std::string filename =
                                "histogram__pop" + boost::lexical_cast<std::string > (pop) +
                                "__" + label_for_dim1_per_pop + boost::lexical_cast<std::string > (jj) +
                                "__" + label_for_dim2_per_pop + boost::lexical_cast<std::string > (kk) +
                                ".dat";
                        _histograms.at(pop, jj, kk) = patterns::createSubject<Histogram_type, sample_type > (sample_space, filename);

                    }

            /* create whole-population histogram matrix */
            {
                for (int ii = 0; ii < dim0_whole; ii++)
                    for (int jj = 0; jj < dim1_whole; jj++) {

                        std::string filename =
                                "histogramWhole__" +
                                label_for_dim1_per_pop +
                                boost::lexical_cast<std::string > (ii) +
                                "__" + label_for_dim2_per_pop + boost::lexical_cast<std::string > (jj) +
                                ".dat";
                        _histograms_whole.at(ii, jj) = patterns::createSubject<Histogram_type, sample_type > (sample_space, filename);
                    }
            }

            /* create histogram observer */
            {
                std::string filename =
                        "histogram__pop" + boost::lexical_cast<std::string > (_pop_to_observe) +
                        "__" + label_for_dim1_per_pop + boost::lexical_cast<std::string > (_jj_per_pop_to_observe) +
                        "__" + label_for_dim2_per_pop + boost::lexical_cast<std::string > (_kk_per_pop_to_observe) +
                        ".log";
                _histogram_observer = patterns::createAttachedObserver<Histogram_Observer_type, Histogram_type > (*_histograms.at(_pop_to_observe, _jj_per_pop_to_observe, _kk_per_pop_to_observe), filename, observer_divisor);
            }

        }

        /** 
         * update sub-population histogram matrices
         */
        void update_histograms(const Pop &pop, const int &jj, const int &kk, const sample_type &sample) {

            _histograms.at(pop.value(), jj, kk)->update(sample);

        }

        /** 
         * update whole-population histogram matrix
         */
        void update_histograms_whole(const int &ii, const int &jj, const sample_type &sample) {

            _histograms_whole.at(ii, jj)->update(sample);

        }


    public:

        /**
         * returns true if statistics of observed sub-population have converged
         */
        virtual const bool converged() const {

            const Histogram_ptr_type observed_histogram = _histograms.at(_pop_to_observe, _jj_per_pop_to_observe, _kk_per_pop_to_observe);
            typedef typename Histogram_type::frequency_t frequency_t;
            const frequency_t threshold_frequency = static_cast<frequency_t> (1.0 / (_error_probability * _error_probability));

            // this assumes that histogram is a CDF or PMF_Bool; 
            // should be generalized so as to use largest_freq if histogram type is PMF or PDF 
            const bool cond1 = observed_histogram->get_end_frequency() >= threshold_frequency;

            typedef typename Histogram_type::number_trials_t number_trials_t;
            const number_trials_t number_trials_max = static_cast<number_trials_t>(1000000000000); // 1e12LL; long long is guaranteed by C++11 to be at least 64 bits
//            const number_trials_t number_trials_max = static_cast<number_trials_t>(1000); // this can be used to end simulation when histogram type is PMF or PDF (see remark about cond1)
            const bool cond2 = observed_histogram->get_number_trials() >= number_trials_max;

            return (cond1 || cond2);

        }

        /** 
         * get histograms for each pop, spe, and node\n
         * \n
         * consolidate methods by returning both sample space and probability distribution in a single object \n
         * this object cannot be a Histogram object because Histogram is noncopyable \n
         * could return ptr to Histogram object but this is bad practice because \n
         * (1) ptr may dangle \n
         * (2) this method could then no longer guarantee not to change internal data as it claims too (const qualifier) \n
         */
        virtual const array::Array3D<SampleSpace_Probability_type> get_results_so_far() const {

            array::Array3D<SampleSpace_Probability_type> array3D_results(_histograms.get_dim0(), _histograms.get_dim1(), _histograms.get_dim2());

            for (int pop = 0; pop < array3D_results.get_dim0(); pop++)
                for (int jj = 0; jj < array3D_results.get_dim1(); jj++)
                    for (int kk = 0; kk < array3D_results.get_dim2(); kk++)
                        array3D_results.at(pop, jj, kk) = fetch_sampleSpace_probability(_histograms.at(pop, jj, kk));

            return array3D_results;
        }

        /** 
         * get whole-population histograms for each spe, and node\n
         */
        virtual const array::Array2D<SampleSpace_Probability_type> get_results_so_far_whole() const {

            array::Array2D<SampleSpace_Probability_type> array2D_results(_histograms_whole.get_dim0(), _histograms_whole.get_dim1());

            for (int ii = 0; ii < array2D_results.get_dim0(); ii++)
                for (int jj = 0; jj < array2D_results.get_dim1(); jj++)
                    array2D_results.at(ii, jj) = fetch_sampleSpace_probability(_histograms_whole.at(ii, jj));

            return array2D_results;
        }


    };



}

#endif	/* DISTRIBUTION_STATISTICS_H */

