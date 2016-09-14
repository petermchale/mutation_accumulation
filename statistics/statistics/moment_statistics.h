#ifndef MOMENT_STATS_H
#define	MOMENT_STATS_H

#include <mutation_accumulation/parameters/parameters_fwd.h> // monte_carlo::Number_Pop

#include "statistics_gatherer.h" // monte_carlo::Statistics_Gatherer

/*************************************************************************/

namespace monte_carlo {

    /**
     * each sub-population is endowed with a matrix of moments\n
     * derive from this to calculate eg MFPT
     */
    template<class Moment_type, class Configuration_type>
    class Moment_Statistics : public Statistics_Gatherer<Configuration_type, typename Moment_type::Results_t> {
    public:

        typedef Moment_type Moment_t;

    private:

        typedef array::Array3D<Moment_type> array3D_moments_type;
        typedef array::Array2D<Moment_type> array2D_moments_type;

        typedef typename Moment_type::sample_t sample_type;
        typedef typename Moment_type::Results_t Results_type;

    private:

        array3D_moments_type _moments; // arbitrary matrix of moments for each sub-population 
        array2D_moments_type _moments_whole; // independent arbitrary matrix of whole-population moments

        const Pop _pop_to_observe;
        const int _jj_per_pop_to_observe;
        const int _kk_per_pop_to_observe;

        const long long int _number_trials_max;

    protected:

        /**
         * constructor \n
         * default ctor of base class is called by compiler
         */
        explicit Moment_Statistics(
                const Number_Pop &number_pop,
                const int &dim1_per_pop,
                const int &dim2_per_pop,
                const int &dim0_whole,
                const int &dim1_whole,
                const Pop &pop_to_observe,
                const int &jj_per_pop_to_observe,
                const int &kk_per_pop_to_observe,
                const long long int &number_trials_max)
        :
        _moments(array3D_moments_type(number_pop.value(), dim1_per_pop, dim2_per_pop)),
        _moments_whole(array2D_moments_type(dim0_whole, dim1_whole)),
        _pop_to_observe(pop_to_observe),
        _jj_per_pop_to_observe(jj_per_pop_to_observe),
        _kk_per_pop_to_observe(kk_per_pop_to_observe),
        _number_trials_max(number_trials_max) {

            /* check template parameter types to supplement "duck typing" */
            BOOST_STATIC_ASSERT((boost::is_base_of<probability::Moment<sample_type>, Moment_type>::value));
        }

        /** 
         * update sub-population histogram matrices
         */
        void update_moments(const Pop &pop, const int &jj, const int &kk, const sample_type &sample) {

            _moments.at(pop.value(), jj, kk).update(sample);

        }

        /** 
         * update whole-population histogram matrix
         */
        void update_moments_whole(const int &ii, const int &jj, const sample_type &sample) {

            _moments_whole.at(ii, jj).update(sample);

        }


    public:

        /**
         * returns true if statistics have converged
         */
        virtual const bool converged() const {

            return _moments.at(_pop_to_observe.value(), _jj_per_pop_to_observe, _kk_per_pop_to_observe).trials() >= _number_trials_max;

        }

        /** 
         * get sub-population probability distributions
         */
        const array::Array3D<Results_type> get_results_so_far() const {

            array::Array3D<Results_type> array3D_results(_moments.get_dim0(), _moments.get_dim1(), _moments.get_dim2());

            for (int pop = 0; pop < array3D_results.get_dim0(); pop++)
                for (int jj = 0; jj < array3D_results.get_dim1(); jj++)
                    for (int kk = 0; kk < array3D_results.get_dim2(); kk++)
                        array3D_results.at(pop, jj, kk) = _moments.at(pop, jj, kk).value();

            return array3D_results;
        }

        /** 
         * get whole-population probability distributions
         */
        const array::Array2D<Results_type> get_results_so_far_whole() const {

            array::Array2D<Results_type> array2D_results(_moments_whole.get_dim0(), _moments_whole.get_dim1());

            for (int ii = 0; ii < array2D_results.get_dim0(); ii++)
                for (int jj = 0; jj < array2D_results.get_dim1(); jj++)
                    array2D_results.at(ii, jj) = _moments_whole.at(ii, jj).value();

            return array2D_results;
        }


    };


}

#endif	/* MOMENT_STATS_H */

