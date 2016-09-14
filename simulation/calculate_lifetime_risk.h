#ifndef CALCULATE_LIFETIME_RISK_H
#define	CALCULATE_LIFETIME_RISK_H

//#define DEBUG_CALCULATE_LIFETIME_RISK_H

#include <cassert> // assert

#include <mutation_accumulation/utility/strings.h> // strings::parse_scalar, etc
#include <mutation_accumulation/configuration/configuration/time_grid.h> // monte_carlo::Uniform_Time_Grid
#include <mutation_accumulation/configuration/configuration/mutation_rates.h> // monte_carlo::MutationRates
#include <mutation_accumulation/configuration/configuration/population2D.h> // monte_carlo::Population2D
#include <mutation_accumulation/statistics/statistics/statistics_lifetime_risk.h> // monte_carlo::Statistics_Lifetime_Risk
#include <mutation_accumulation/configuration/configuration/path_policy.h> // monte_carlo::Wait_For_Last_Species_In_All_SubPops
#include <mutation_accumulation/simulation/generate_statistics.h> // monte_carlo::generate_statistics
#include <mutation_accumulation/utility/grid.h> // grid::make_logarithmic_grid
#include <mutation_accumulation/configuration/utilities/create.h> // monte_carlo::create_configuration


/*************************************************************************/

namespace monte_carlo {

    namespace calculate_lifetime_risk_detail {

        typedef double error_type;
        typedef int divisor_type;

        /** 
         * calculate lifetime risk for a particular parameter set and a generic homeostatic stochastic process 
         * assumes that number of species = 1 + (number of mutation rates), which isn't true for "diamond"; diamond configuration will complain at compile time
         */
        template <class Configuration_Policy>
        const double calculate_lifetime_risk(
                const typename Configuration_Policy::population_t &N0,
                const MutationRates &uu,
                const Symmetry &symmetry,
                const typename Configuration_Policy::time_t &time_span_path,
                const error_type &error_probability,
                const divisor_type &observer_divisor) {

            typedef typename Configuration_Policy::time_t time_type;
            typedef typename Configuration_Policy::population_t population_type;

            /* create a uniform grid of time points at which to sample configuration */
            Uniform_Time_Grid<time_type> time_grid(time_span_path);

            /* create initial population matrix with a size determined by size of uu */
            const int number_species = uu.size() + 1;
            std::vector<population_type> NN_vec(number_species);
            NN_vec.at(0) = N0;
            for (int ii = 1; ii < NN_vec.size(); ii++)
                NN_vec.at(ii) = static_cast<population_type> (0);
            const Population2D<population_type> NN(NN_vec);
#ifdef DEBUG_CALCULATE_LIFETIME_RISK_H
            /* check NN */
            std::cout << "NN = ";
            for (int pop = 0; pop < NN.number_sub_pops(); pop++) {
                for (int spe = 0; spe < NN.number_species(); spe++)
                    std::cout << "<" << NN.at(pop, spe) << ">" << " ";
                std::cout << std::endl;
            }
#endif

#ifdef DEBUG_CALCULATE_LIFETIME_RISK_H
            /* compare uu with NN */
            assert(uu.get_dim() == (NN.number_species() - 1));
#endif

            /* initialize configuration */
            typename Configuration_Policy::category configuration_category;
            const Configuration_Policy configuration_init = create_configuration<Configuration_Policy > (configuration_category, NN, uu, symmetry, time_grid);

            /* observe statistics in first sub-population */
            const int first_pop = 0;
            const Pop pop_to_observe(first_pop);

            /* observe statistics in last species */
            const int last_species = NN.number_species() - 1;
            const Spe spe_to_observe(last_species);

            /* create statistics gatherer type */
            typedef Statistics_Lifetime_Risk<Configuration_Policy> Statistics_Policy;

            /* instantiate statistics gatherer */
            Statistics_Policy statistics(
                    Number_Pop(NN.number_sub_pops()),
                    Number_Spe(NN.number_species()),
                    error_probability,
                    pop_to_observe,
                    spe_to_observe,
                    observer_divisor);

            /* path policy */
            typedef Wait_For_Last_Species_In_All_SubPops<Configuration_Policy> Path_Policy;

            /* do Monte Carlo simulation; gather statistics; find successful paths */
            Generate_Statistics<Path_Policy, Configuration_Policy, Statistics_Policy >::implement(configuration_init, statistics);

            /* calculate lifetime risk and return */
            return probability_mutation_fate(statistics);

        }

        /** 
         * loop over population sizes (and symmetry values) and print out lifetime risk using generic homeostatic stochastic process
         */
        template <class Configuration_Policy>
        const void calculate_lifetime_risk__loop_over_N() {

            typedef typename Configuration_Policy::time_t time_type;
            typedef typename Configuration_Policy::population_t population_type;

            /* open input file */
            const std::string filename = "main.in";
            boost::shared_ptr<std::ifstream> ifstream_ptr = monte_carlo::open_file_for_input(filename);

            /* read in N parameters */
            const int number_N = strings::parse_scalar<int > (*ifstream_ptr);
            const population_type N_lower = strings::parse_scalar<population_type > (*ifstream_ptr);
            const population_type N_upper = strings::parse_scalar<population_type > (*ifstream_ptr);

            /* read in mutation rates */
            const MutationRates uu(strings::parse_vector<MutationRates::data_t > (*ifstream_ptr));
#ifdef DEBUG_CALCULATE_LIFETIME_RISK_H 
            /* check uu */
            std::cout << "uu = ";
            for (int ii = 0; ii < uu.get_dim(); ii++)
                std::cout << "<" << uu.at(ii) << "> ";
            std::cout << std::endl;
#endif

            /* read symmetry strings */
            const std::vector<std::string> symmetry_strings = strings::read_vectorStrings(*ifstream_ptr, "&");
#ifdef DEBUG_CALCULATE_LIFETIME_RISK_H 
            /* check ss */
            std::cout << "ss strings = ";
            for (int ii = 0; ii < symmetry_strings.size(); ii++)
                std::cout << "<" << symmetry_strings.at(ii) << "> ";
            std::cout << std::endl;
#endif

            /* convert symmetry strings to values */
            const std::vector<Symmetry> symmetry_values = strings::vectorStrings_to_vectorValues<Symmetry > (symmetry_strings);
#ifdef DEBUG_CALCULATE_LIFETIME_RISK_H 
            /* check ss */
            std::cout << "ss values = ";
            for (int ii = 0; ii < symmetry_values.size(); ii++)
                std::cout << "<" << symmetry_values.at(ii).value() << "> ";
            std::cout << std::endl;
#endif

            /* read in time span for path */
            const time_type time_span_path = strings::parse_scalar<time_type > (*ifstream_ptr);

            /* read in error in cumulative probability at population span */
            const error_type error_probability = strings::parse_scalar<error_type > (*ifstream_ptr);

            /* write log data every observer_divisor times observer is notified */
            const divisor_type observer_divisor = strings::parse_scalar<divisor_type > (*ifstream_ptr);

            {
                /* create file to store symmetry strings for matlab use */
                boost::shared_ptr<std::ofstream> ofstream_ptr_ssFile = monte_carlo::open_file_for_output("ss.in");

                /* create list of logarithmically spaced population sizes */
                const std::vector<population_type> N_grid = grid::make_logarithmic_grid(number_N, N_lower, N_upper);

                for (int ii = 0; ii < symmetry_values.size(); ii++) {

                    /* record symmetry string */
                    *ofstream_ptr_ssFile << symmetry_strings.at(ii) << std::endl;

                    /* create file to store lifetime risks for various N */
                    const std::string filename = "lifetime_risk__s" + symmetry_strings.at(ii) + ".dat";
                    boost::shared_ptr<std::ofstream> ofstream_ptr = monte_carlo::open_file_for_output(filename);

                    for (int jj = 0; jj < N_grid.size(); jj++) {

                        *ofstream_ptr << std::setw(30) << std::setprecision(20) << N_grid.at(jj);
                        *ofstream_ptr << std::setw(30) << std::setprecision(20) << calculate_lifetime_risk<Configuration_Policy > (N_grid.at(jj), uu, symmetry_values.at(ii), time_span_path, error_probability, observer_divisor);
                        *ofstream_ptr << std::endl;
                    }
                }

                /* indicate that simulation has finished */
                monte_carlo::done();

            }
        }

        /** 
         * read in a particular parameter set; calculate lifetime risk; dump result to disk 
         */
        template <class Configuration_Policy>
        const void lifetime_risk__read_calculate_dump() {

            typedef typename Configuration_Policy::time_t time_type;
            typedef typename Configuration_Policy::population_t population_type;

            /* open input file */
            const std::string filename = "main.in";
            boost::shared_ptr<std::ifstream> ifstream_ptr = monte_carlo::open_file_for_input(filename);

            /* read in initial population sizes */
            typedef Population2D<population_type> Population2D_type;
            const Population2D_type NN(strings::parse_matrix<typename Population2D_type::population_t > (*ifstream_ptr));

            /* read in time span */
            const time_type time_span = strings::parse_scalar<time_type > (*ifstream_ptr);

            /* read in mutation rates */
            const MutationRates uu(strings::parse_vector<MutationRates::data_t > (*ifstream_ptr));

            /* compare uu with NN */
            assert(uu.get_dim() == (NN.number_species() - 1));

            /* read in symmetry */
            const Symmetry symmetry = strings::parse_scalar<Symmetry > (*ifstream_ptr);

            /* read in error in end of life risk */
            const error_type error_probability = strings::parse_scalar<error_type > (*ifstream_ptr);

            /* write log data every observer_divisor times observer is notified */
            const divisor_type observer_divisor = strings::parse_scalar<divisor_type > (*ifstream_ptr);

            {

                /* calculate and dump lifetime risk to disk */
                boost::shared_ptr<std::ofstream> ofstream_ptr = monte_carlo::open_file_for_output("lifetime_risk.dat");
                *ofstream_ptr << std::setw(30) << std::setprecision(20);
                *ofstream_ptr << calculate_lifetime_risk<Configuration_Policy > (NN.at(0, 0), uu, symmetry, time_span, error_probability, observer_divisor);
                *ofstream_ptr << std::endl;

                /* indicate that simulation has finished */
                monte_carlo::done();

            }
        }
    }

    using calculate_lifetime_risk_detail::calculate_lifetime_risk__loop_over_N;
    using calculate_lifetime_risk_detail::lifetime_risk__read_calculate_dump;

}



#endif	/* CALCULATE_LIFETIME_RISK_H */

