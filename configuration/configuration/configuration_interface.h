#ifndef CONFIG_INTERFACE_H
#define	CONFIG_INTERFACE_H

#include <mutation_accumulation/parameters/parameters_fwd.h> // monte_carlo::Pop, etc

#include "random_fwd.h" // monte_carlo::base_generator_type

/*************************************************************************/

namespace monte_carlo {

    /**
     * pure abstract base class that defines the interface of Configuration\n
     * may use this class polymorphically (virtual destructor) \n
     * may not instantiate (ctor is protected)\n
     */
    template <class time_type, class population_type>
    class Configuration_Interface {
    public:

        typedef time_type time_t;
        typedef population_type population_t;

    protected:

        /**
         * custom constructor
         * this is an empty constructor\n
         * compiler should generate code to construct private data members (Item 30)
         */
        explicit Configuration_Interface() {

        }


    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Configuration_Interface() {

        }

        /**
         * get current time
         */
        virtual const time_type get_time() const = 0; 

        /**
         * get number of type-spe cells in sub-population pop
         */
        virtual const population_type get_population(const Pop &pop, const Spe &spe) const = 0; 

        /**
         * number of sub-populations
         */
        virtual const int number_sub_pops() const = 0;

        /**
         * number of species
         */
        virtual const int number_species() const = 0;

        /**
         * number of filled nodes
         */
        virtual const int number_filled_nodes() const = 0;

        /**
         * get time at which sub-population pop extinguished
         */
        virtual const time_type extinction_times(const Pop & pop) const = 0;

        /**
         * get time at which first spe-type cell appears in sub-population pop
         */
        virtual const time_type mutation_times(const Pop &pop, const Spe & spe) const = 0;

        /**
         * get path time \n
         */
        virtual const time_type get_path_time(const Node &node) const = 0;

        /**
         * get time of last node in path\n
         */
        virtual const time_type get_last_node_time() const = 0;

        /**
         * get path population \n
         */
        virtual const population_type get_path_population(const Pop &pop, const Spe &spe, const Node &node) const = 0;

        /**
         * transition to next state of stochastic process
         */
        virtual void transition(base_generator_type &base_rand_gen) = 0; 


        /**
         * determine if end of trajectory has been reached
         */
        virtual const bool end() const = 0;


    };



}

#endif	

