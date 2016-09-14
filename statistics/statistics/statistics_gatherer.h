#ifndef STATISTICS_GATHERER_H
#define	STATISTICS_GATHERER_H

#include <boost/type_traits.hpp> // boost::is_base_of

#include <mutation_accumulation/configuration/configuration/configuration_interface.h> // monte_carlo::Configuration_Interface 
#include <mutation_accumulation/array/array3D.h> // array::Array3D 
#include <mutation_accumulation/array/array2D.h> // array::Array2D 

/*************************************************************************/

namespace monte_carlo {

    /**
     * abstract base class that stores statistics generated by Configuration class\n
     * defines the interface used by generate_statistics(..) \n
     * may use this class polymorphically (virtual destructor) \n
     * may not instantiate (ctor is protected)\n
     */
    template <class Configuration_type, class Results_type>
    class Statistics_Gatherer {
    public:

        typedef Configuration_type Configuration_t;
        typedef Results_type Results_t;
        
    protected:

        explicit Statistics_Gatherer() {

            /* check template parameter types to supplement "duck typing" */
            typedef typename Configuration_type::time_t time_type; // Item 42
            typedef typename Configuration_type::population_t population_type;
            BOOST_STATIC_ASSERT((boost::is_base_of<Configuration_Interface<time_type, population_type>, Configuration_type>::value));
        }

    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Statistics_Gatherer() {

        }


        /** 
         * dump results of a particular trial\n
         * cannot pull the loop over pop, spe, node up to this level without also pulling up the size of the 3D array\n
         * (statistics for extinction times makes the problem concrete)
         */
        virtual void dump(const Configuration_type &configuration) = 0;

        /**
         * returns true if statistics of observed population have converged
         */
        virtual const bool converged() const = 0;

        /** 
         * return results as 3D array so that number of pops, spe, and nodes can be retrieved \n
         */
        virtual const array::Array3D<Results_type> get_results_so_far() const = 0;

        /** 
         * return whole-population results as 2D array so that number of spe, and nodes can be retrieved \n
         */
        virtual const array::Array2D<Results_type> get_results_so_far_whole() const = 0;



    };



}

#endif	/* STATISTICS_GATHERER_H */

