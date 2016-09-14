#ifndef PRINT_TRAJECTORIES_H
#define	PRINT_TRAJECTORIES_H

#include <boost/lexical_cast.hpp> //boost::lexical_cast
#include <boost/static_assert.hpp> // BOOST_STATIC_ASSERT
#include <boost/type_traits/is_base_of.hpp> //  boost::is_base_of

#include <mutation_accumulation/configuration/configuration/random_fwd.h> // monte_carlo::base_generator_type
#include <mutation_accumulation/configuration/configuration/configuration_interface.h> // monte_carlo::Configuration_Interface
#include <mutation_accumulation/configuration/configuration/path_policy.h> // monte_carlo::Path_Policy_Base
#include <mutation_accumulation/simulation/files.h> // monte_carlo::open_file_for_output 
#include <mutation_accumulation/configuration/configuration/branching_discrete.h> // monte_carlo::Branching_Discrete
#include <mutation_accumulation/configuration/configuration/moran4.h> // monte_carlo::Moran4
#include <mutation_accumulation/configuration/utilities/print.h> // <<configuration; print_path

#define PRINT_TRANSITIONS // don't know how to make PRINT_TRANSITIONS a run-time flag specific to Moran4 (vs Branching_Discrete)

/*************************************************************************/

namespace monte_carlo {

    namespace print_trajectories_detail {

        /** 
         * check template parameters of print functions
         */
        template <class Configuration_Policy, class Path_Policy>
        void check_template_parameters() {

            /* typedefs */
            typedef typename Configuration_Policy::time_t time_t;
            typedef typename Configuration_Policy::population_t population_t;

            /* check policy type to supplement "duck typing" */
            BOOST_STATIC_ASSERT((boost::is_base_of<Configuration_Interface<time_t, population_t>, Configuration_Policy>::value));
            BOOST_STATIC_ASSERT((boost::is_base_of<Path_Policy_Base<Configuration_Policy>, Path_Policy>::value));

        }

        /** 
         * print single trajectory at nodes and transitions
         */
        template <class Configuration_Policy, class Path_Policy>
        void print_trajectory(
        const Configuration_Policy &configuration_init,
        const Path_Policy &path_policy,
        const std::string &file_name,
        base_generator_type &base_rand_gen) {

#ifdef PRINT_TRANSITIONS
            /* open file to store trajectory at transitions */
            boost::shared_ptr<std::ofstream> ofstream_transitions_ptr = open_file_for_output(file_name + "_transitions");
#endif

            /* initialize the state of the stochastic process with a deterministic configuration\n
             * could initialize with a random configuration */
            Configuration_Policy configuration(configuration_init);
#ifdef PRINT_TRANSITIONS
            *ofstream_transitions_ptr << configuration << std::endl;
#endif

            /* generate trajectory of the stochastic process */
            while (!path_policy.terminate(configuration)) {
                configuration.transition(base_rand_gen);
#ifdef PRINT_TRANSITIONS
                *ofstream_transitions_ptr << configuration << std::endl;
#endif
            }

            /* print trajectory at nodes */
            print_path(configuration, file_name + "_nodes");

        }
    }

    /**
     * print successful trajectories (defined by path policy) using stochastic process (defined by configuration policy)\n
     * strategy/policy design pattern\n
     * static (compile-time) polymorphism via templates
     */
    template <class Configuration_Policy, class Path_Policy>
    void print_successful_trajectories(const Configuration_Policy &configuration_init, const Path_Policy &path_policy) {

        /* check template parameters */
        print_trajectories_detail::check_template_parameters<Configuration_Policy, Path_Policy > ();

        /* perform small number of a-posteriori successful trials */
        for (int trial = 0; trial < path_policy.successful_PRNG_states().size(); trial++) {

            /* copy construct random number engine */
            base_generator_type base_rand_gen(path_policy.successful_PRNG_states().at(trial)); // std::time(0) is a random seed

            std::string file_name = "successful_trajectory" + boost::lexical_cast<std::string > (trial) + ".dat";
            print_trajectories_detail::print_trajectory(configuration_init, path_policy, file_name, base_rand_gen);

        }

    }

    /**
     * print random trajectories using stochastic process (defined by configuration policy)\n
     * path policy is used only to terminate path\n
     * strategy/policy design pattern\n
     * static (compile-time) polymorphism via templates
     */
    template <class Configuration_Policy, class Path_Policy>
    void print_random_trajectories(const Configuration_Policy &configuration_init, const Path_Policy &path_policy) {

        /* check template parameters */
        print_trajectories_detail::check_template_parameters<Configuration_Policy, Path_Policy > ();

        /* random number engine */
        base_generator_type base_rand_gen(std::time(0)); // static_cast<unsigned int> (0)

        /* perform small number of random trials */
        for (int trial = 0; trial < 10; trial++) {

            std::string file_name = "random_trajectory" + boost::lexical_cast<std::string > (trial) + ".dat";
            print_trajectories_detail::print_trajectory(configuration_init, path_policy, file_name, base_rand_gen);

        }

    }


}



#endif	/* PRINT_TRAJECTORIES_H */

