/* 
 * Author: petermchale
 *
 * Created on June 18, 2012, 3:57 PM
 */

#ifndef DERIVED_CLASSES_H
#define	DERIVED_CLASSES_H

#include <mutation_accumulation/utility/data_traits.h> // data_types::continuous_type

#include <limits> // std::numeric_limits
#include <cassert> // assert

#include "parameter_base.h"

/*************************************************************************/

/* consider moving implementation of implicitly inlined constructors to implementation file
 * (Item 30) */

namespace parameters {

    namespace ranges {
        
        /* allowable parameter ranges */
        const double symmetry_min = 0.0;
        const double symmetry_max = 1.0;

        const double symmetric_renewal_min = 0.0;
        const double symmetric_renewal_max = 1.0;

        const int number_pop_min = 0;
        const int number_pop_max = 100;

        const int change_pop_size_min = std::numeric_limits<int>::min();
        const int change_pop_size_max = std::numeric_limits<int>::max();

        const int number_spe_min = 0;
        const int number_spe_max = 100;

        const int number_sampled_times_min = 0; 
        const int number_sampled_times_max = static_cast<int>(1e5); 
        
        /**
         * check that model parameter value lies in allowable range
         */
        inline void check_value(const int &value__, const int &value_min__, const int &value_max__) {

            assert(value__ >= value_min__);
            assert(value__ <= value_max__);

        }

        /**
         * check that model parameter value lies in allowable range 
         */
        inline void check_value(const double &value__, const double &value_min__, const double &value_max__) {

            assert(value__ > value_min__ - 1e-8);
            assert(value__ < value_max__ + 1e-8);

        }

    }

    /**
     * symmetry parameter
     */
    class Symmetry : public Parameter_base<double> {
    public: 
        typedef data_types::continuous_type category;

    public:

        /**
         * ctor
         */
        explicit Symmetry(const double &value__) : Parameter_base<double>(value__, ranges::symmetry_min, ranges::symmetry_max) {

            ranges::check_value(value_, value_min_, value_max_);
            
        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit Symmetry() : Parameter_base<double>(0.0, ranges::symmetry_min, ranges::symmetry_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }


    };


    /**
     * symmetry parameter
     */
    class SymmetricRenewal : public Parameter_base<double> {
    public:

        /**
         * ctor
         */
        explicit SymmetricRenewal(const double &value__) : Parameter_base<double>(value__, ranges::symmetric_renewal_min, ranges::symmetric_renewal_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit SymmetricRenewal() : Parameter_base<double>(0.5, ranges::symmetric_renewal_min, ranges::symmetric_renewal_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }


    };



    /**
     * number of sub-populations parameter
     */
    class Number_Pop : public Parameter_base<int> {
    public:

        /**
         * ctor
         */
        explicit Number_Pop(const int &value__) : Parameter_base<int>(value__, ranges::number_pop_min, ranges::number_pop_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit Number_Pop() : Parameter_base<int>(0, ranges::number_pop_min, ranges::number_pop_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

    };

    /**
     * population index
     */
    class Pop : public Parameter_base<int> {
    public:

        /**
         * ctor
         */
        explicit Pop(const int &value__) : Parameter_base<int>(value__, ranges::number_pop_min, ranges::number_pop_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

    };

    /**
     * number of cells to add (positive value__) or take away (negative value__) from a population
     */
    class Change_Pop_Size : public Parameter_base<int> {
    public:

        /**
         * ctor
         */
        explicit Change_Pop_Size(const int &value__) : Parameter_base<int>(value__, ranges::change_pop_size_min, ranges::change_pop_size_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

    };


    /**
     * number of species parameter
     */
    class Number_Spe : public Parameter_base<int> {
    public:

        /**
         * ctor
         */
        explicit Number_Spe(const int &value__) : Parameter_base<int>(value__, ranges::number_spe_min, ranges::number_spe_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

    };

    /**
     * species index
     */
    class Spe : public Parameter_base<int> {
    public:

        /**
         * ctor
         */
        explicit Spe(const int &value__) : Parameter_base<int>(value__, ranges::number_spe_min, ranges::number_spe_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

    };

    
    /**
     * number of time points at which to sample a path
     */
    class Number_Nodes : public Parameter_base<int> {
    public:

        /**
         * ctor
         */
        explicit Number_Nodes(const int &value__) : Parameter_base<int>(value__, ranges::number_sampled_times_min, ranges::number_sampled_times_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

    };

    /**
     * time index for sampling a path
     */
    class Node : public Parameter_base<int> {
    public:

        /**
         * ctor
         */
        explicit Node(const int &value__) : Parameter_base<int>(value__, ranges::number_sampled_times_min, ranges::number_sampled_times_max) {

            ranges::check_value(value_, value_min_, value_max_);

        }

    };


} 

#endif	/* DERIVED_CLASSES_H */

