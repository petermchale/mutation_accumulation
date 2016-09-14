#ifndef CALCULATE_HISTOGRAM_TRAJS_H
#define	CALCULATE_HISTOGRAM_TRAJS_H

#include <mutation_accumulation/statistics/statistics/statistics_mutation.h> // monte_carlo::Statistics_Mutation
#include <mutation_accumulation/configuration/configuration/path_policy.h> // monte_carlo::Wait_For_Last_Species_In_All_SubPops
#include <mutation_accumulation/simulation/generate_statistics.h> // monte_carlo::generate_statistics
#include <mutation_accumulation/simulation/print_trajectories.h> // monte_carlo::print_xxx_trajectories
#include <mutation_accumulation/configuration/utilities/create_configuration.h> // monte_carlo::create_configuration
#include <mutation_accumulation/simulation/read_policy.h>  
#include <mutation_accumulation/statistics/utilities/create_statistics.h>

/*************************************************************************/

namespace monte_carlo {

    /** 
     * calculate histogram and sample trajectories for an arbitrary stochastic process (not necessarily homeostatic)
     * Raw_Data_Policy defines whether and how to print out random samples \n
     * function template not used because default template parameters not allowed in function templates
     * If the default is specified for a template parameter, each subsequent template parameter must have a default argument
     */
    template <
    class Histogram_Policy,
    class Configuration_Policy,
    template <class Configuration_type> class Raw_Data_Policy = Raw_Data_Null,
    template <class Configuration_type> class Read_Policy = Read_Homeostasis_Policy    
    >
    class Calculate_Histogram_Trajs {
    public:

        static const void implement() {

            /* set up policy object that reads input */
            Read_Policy<Configuration_Policy> read_policy("main.in");

            {
                /* initialize configuration */
                typename Configuration_Policy::category configuration_category;
                const Configuration_Policy configuration_init = create_configuration<Configuration_Policy > (configuration_category, read_policy);

                /* create statistics gatherer type */
                typedef Statistics_Mutation<Histogram_Policy, Configuration_Policy> Statistics_Policy;

                /* create instance of a statistics gatherer */
                Statistics_Policy statistics = create_statistics<Statistics_Policy, Configuration_Policy > (read_policy);

                /* path policy */
                typedef Wait_For_Last_Species_In_All_SubPops<Configuration_Policy> Path_Policy;

                /* do Monte Carlo simulation; gather statistics; find successful paths */
                Path_Policy path_policy = Generate_Statistics<Path_Policy, Configuration_Policy, Statistics_Policy, Raw_Data_Policy >::implement(configuration_init, statistics);

                /* print successful trajectories */
                print_successful_trajectories(configuration_init, path_policy);

                /* print random trajectories */
                print_random_trajectories(configuration_init, path_policy);

                /* indicate that simulation has finished */
                done();

            }
        }
    };


}



#endif	/* CALCULATE_HISTOGRAM_TRAJS_H */

