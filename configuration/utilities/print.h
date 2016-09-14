#ifndef PRINT_CONFIGURATION_H
#define	PRINT_CONFIGURATION_H

#include <iomanip> // std::setw

#include <mutation_accumulation/configuration/configuration/configuration_interface.h> // monte_carlo::Configuration_Interface
#include <mutation_accumulation/simulation/files.h> // monte_carlo::open_file_for_output 

/*************************************************************************/

namespace monte_carlo {


    namespace print_utility {
        
        /* if one of these two overloaded functions is deleted the other will take over by up- or down-casting tt */

        inline void print_time(std::ostream &out_stream, const int &tt) {

            out_stream << std::setw(10) << tt;

        }

        inline void print_time(std::ostream &out_stream, const double &tt) {

            out_stream << std::setw(20) << std::setprecision(10) << tt;

        }

    }

    /**
     * print out state of Configuration_Interface object
     */
    template <class time_type, class population_type>
    std::ostream & operator<<(std::ostream &out_stream, const Configuration_Interface<time_type, population_type> &configuration) {

        print_utility::print_time(out_stream, configuration.get_time());

        for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
            for (int spe = 0; spe < configuration.number_species(); spe++)
                out_stream << std::setw(20) << configuration.get_population(Pop(pop), Spe(spe)) << " ";

        return out_stream;

    }

    /**
     * print out time course of path object in Configuration_Interface object
     */
    template <class time_type, class population_type>
    void print_path(const Configuration_Interface<time_type, population_type> &configuration, const std::string &file_name) {

        boost::shared_ptr<std::ofstream> ofstream_ptr = open_file_for_output(file_name);

        for (int node = 0; node < configuration.number_filled_nodes(); node++) {
            print_utility::print_time(*ofstream_ptr, configuration.get_path_time(Node(node)));

            for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
                for (int spe = 0; spe < configuration.number_species(); spe++)
                    *ofstream_ptr << std::setw(20) << configuration.get_path_population(Pop(pop), Spe(spe), Node(node)) << " ";
            
            *ofstream_ptr << std::endl;
        }
    }


}

#endif	/* PRINT_CONFIGURATION_H */

