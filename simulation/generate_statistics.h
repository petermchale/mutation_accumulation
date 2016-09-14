#ifndef GENERATE_STATISTICS_H
#define	GENERATE_STATISTICS_H

//#define DEBUG_GENERATE_STATISTICS

#include <boost/static_assert.hpp> // BOOST_STATIC_ASSERT
#include <boost/type_traits/is_base_of.hpp> //  boost::is_base_of

#include <mutation_accumulation/configuration/configuration/random_fwd.h> // monte_carlo::base_generator_type
#include <mutation_accumulation/configuration/configuration/configuration_interface.h> // monte_carlo::Configuration_Interface
#include <mutation_accumulation/statistics/statistics/statistics_gatherer.h> // monte_carlo::Statistics_Gatherer

#include "raw_data.h"

/*************************************************************************/

namespace monte_carlo {

    /** 
     *  print debug info
     */
    template <class Configuration_Policy>
    void print_debug_info(const Configuration_Policy &configuration) {

        std::cout << "***************************************" << std::endl;

        std::cout << "time = " << configuration.get_time() << std::endl;
        std::cout << std::endl;

        std::cout << "configuration = " << std::endl;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++) {
            for (int spe = 0; spe < configuration.number_species(); spe++)
                std::cout << configuration.get_population(Pop(pop), Spe(spe)) << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;

        int total_population_size = 0;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++) {
            for (int spe = 0; spe < configuration.number_species(); spe++)
                total_population_size += configuration.get_population(Pop(pop), Spe(spe));
        }
        std::cout << "total population size = " << total_population_size << std::endl;
        std::cout << std::endl;

        std::cout << "configuration::mutation_times = " << std::endl;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++) {
            for (int spe = 0; spe < configuration.number_species(); spe++)
                std::cout << configuration.mutation_times(Pop(pop), Spe(spe)) << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "configuration::extinction_times = " << std::endl;
        for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
            std::cout << configuration.extinction_times(Pop(pop)) << std::endl;
        std::cout << extinction_time_whole(configuration) << std::endl;
        std::cout << std::endl;

        std::cout << "***************************************" << std::endl;


    }


    /** 
     * calculate statistics (defined by statistics policy) of stochastic process (defined by configuration policy)\n
     * path policy defines whether given stochastic path is successful\n
     * Raw_Data_Policy defines whether and how to print out random samples \n
     * strategy/policy design pattern\n
     * static (compile-time) polymorphism via templates \n
     * function template not used because default template parameters not allowed in function templates
     */
    template <class Path_Policy, class Configuration_Policy, class Statistics_Policy, template <class Configuration_type> class Raw_Data_Policy = Raw_Data_Null>
    class Generate_Statistics {
    public:

        static const Path_Policy implement(const Configuration_Policy &configuration_init, Statistics_Policy &statistics) {

            /* typedefs */
            typedef typename Configuration_Policy::time_t time_t;
            typedef typename Configuration_Policy::population_t population_t;
            typedef typename Statistics_Policy::Results_t Results_t;

            /* check policy type to supplement "duck typing" */
            BOOST_STATIC_ASSERT((boost::is_base_of<Path_Policy_Base<Configuration_Policy>, Path_Policy>::value));
            BOOST_STATIC_ASSERT((boost::is_base_of<Configuration_Interface<time_t, population_t>, Configuration_Policy>::value));
            BOOST_STATIC_ASSERT((boost::is_base_of<Statistics_Gatherer<Configuration_Policy, Results_t>, Statistics_Policy>::value));
            BOOST_STATIC_ASSERT((boost::is_base_of<Raw_Data_Policy_Base<Configuration_Policy>, Raw_Data_Policy<Configuration_Policy> >::value));

            /* PRNG */
            base_generator_type base_rand_gen(std::time(0)); // static_cast<unsigned int> (0)

            /* path policy object determines how to handle the paths */
            Path_Policy path_policy;

            /* set up policy object that handles the printing of raw data */
            Raw_Data_Policy<Configuration_Policy> raw_data_policy;

            /* perform multiple trials of the stochastic process */
            while (!statistics.converged()) {

                /* record initial state of random number generator */
                base_generator_type base_rand_gen_init(base_rand_gen);

                /* initialize the state of the stochastic process with a deterministic configuration\n
                 * could initialize with a random configuration */
                Configuration_Policy configuration(configuration_init);

                /* generate trajectory of the stochastic process */
                while (!path_policy.terminate(configuration))
                    configuration.transition(base_rand_gen);

#ifdef DEBUG_GENERATE_STATISTICS
                print_debug_info(configuration);
#endif

                /* gather statistics from the stochastic process */
                statistics.dump(configuration);

                /* print raw data to disk */
                raw_data_policy.print(configuration);

                /* save initial PRNG state if path was successful */
                path_policy.record_success(configuration, base_rand_gen_init);

            }

            return path_policy;

        }

    };


}



#endif	/* GENERATE_STATISTICS_H */

