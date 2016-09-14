#ifndef CONFIG_BASE_H
#define	CONFIG_BASE_H

#include <mutation_accumulation/utility/data_traits.h> // data_types::discrete_type, etc

#include "configuration_interface.h" // monte_carlo::Configuration_Interface
#include "extinction_times.h" // monte_carlo::ExtinctionTimes
#include "mutation_times.h" // monte_carlo::MutationTimes
#include "path.h" // monte_carlo::Path

/*************************************************************************/

namespace monte_carlo {

    /**
     * abstract base class that stores the configuration of arbitrary number of sub-populations containing arbitrary number of species\n
     * \n
     * may use this class polymorphically (virtual destructor) \n
     * may not instantiate (ctor is protected)\n
     */
    template <class time_type, class population_type>
    class Configuration : public Configuration_Interface<time_type, population_type> {
    private:

        /* state variables required to implement non-member non-friend print function */
        time_type _time;

        typedef Population2D<population_type> population2D_type;
        population2D_type _population2D;
        population2D_type _population2D_old;

        /* variables required to implement required path-dependent methods: 
         * extinction_times(..), mutation_times(..), path(..) */
        ExtinctionTimes<time_type, population_type> _extinction_times;

        MutationTimes<time_type, population_type> _mutation_times;

        typedef Path<time_type, population_type> path_t;
        path_t _path;


    private:

        /**
         * update populations and time
         */
        virtual void update_populations_and_time(base_generator_type &base_rand_gen) = 0;

        /**
         * update path; discrete time
         */
        void update_path(data_types::discrete_type) {

            _path.update(_time, _population2D);
        }

        /**
         * update path; continuous time
         */
        void update_path(data_types::continuous_type) {

            _path.update(_time, _population2D, _population2D_old);
        }


    protected:

        /**
         * custom constructor
         * interface-class default ctor called by compiler
         */
        explicit Configuration(
                const population2D_type &population2D,
                const Uniform_Time_Grid<time_type> &time_grid)
        :
        _time(static_cast<time_type> (0)),
        _population2D(population2D),
        _population2D_old(population2D),
        _extinction_times(population2D.number_sub_pops()),
        _mutation_times(population2D),
        _path(path_t(population2D, time_grid)) {

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist\n
         * interface-class default ctor called by compiler
         * this is an empty constructor\n
         * compiler should generate code to construct private data members (Item 30)
         */
        explicit Configuration() {

        }

        /**
         * set current time
         */
        void set_time(const time_type &new_time) {

            _time = new_time;
        }

        /**
         * set number of type-spe cells in sub-population pop
         */
        void set_population(const Pop &pop, const Spe &spe, const population_type &new_population_size) {

            _population2D.at(pop.value(), spe.value()) = new_population_size;
        }


    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Configuration() {

        }

        /**
         * get current time
         */
        virtual const time_type get_time() const {

            return _time;
        }

        /**
         * get number of type-spe cells in sub-population pop
         */
        virtual const population_type get_population(const Pop &pop, const Spe &spe) const {

            return _population2D.at(pop.value(), spe.value());
        }

        /**
         * number of sub-populations
         */
        virtual const int number_sub_pops() const {

            return _population2D.number_sub_pops();
        }

        /**
         * number of species
         */
        virtual const int number_species() const {

            return _population2D.number_species();
        }

        /**
         * number of filled nodes
         */
        virtual const int number_filled_nodes() const {

            return _path.number_filled_nodes();
        }

        /**
         * get time at which sub-population pop extinguished
         */
        virtual const time_type extinction_times(const Pop & pop) const {

            return _extinction_times(pop);
        }

        /**
         * get time at which first spe-type cell appears in sub-population pop
         */
        virtual const time_type mutation_times(const Pop &pop, const Spe & spe) const {

            return _mutation_times(pop, spe);
        }

        /**
         * get path time \n
         */
        virtual const time_type get_path_time(const Node &node) const {
            
            return _path.get_time(node);
        }

        /**
         * get time of last node in path\n
         */
        virtual const time_type get_last_node_time() const {
            
            return _path.last_node_time();
        }

        /**
         * get path population \n
         */
        virtual const population_type get_path_population(const Pop &pop, const Spe &spe, const Node &node) const {
            
            return _path.at(pop, spe, node);
        }

        /**
         * transition to next state of stochastic process
         */
        virtual void transition(base_generator_type &base_rand_gen) {

            /* store old value of population matrix before it is changed */
            _population2D_old = _population2D;

            /* update populations and time */
            update_populations_and_time(base_rand_gen);

            /* update extinction times */
            _extinction_times.update(_time, _population2D);

            /* update mutation times */
            _mutation_times.update(_time, _population2D);

            /* update path */
            typename data_types::data_traits<time_type>::category time_category;
            update_path(time_category);

        }

        /**
         * determine if end of trajectory has been reached
         */
        virtual const bool end() const {

            return _path.complete();

        }



    };



}

#endif	

