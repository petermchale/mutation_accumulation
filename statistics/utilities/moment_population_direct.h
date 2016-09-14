#ifndef MOMENT_POPULATION_DIRECT_H
#define	MOMENT_POPULATION_DIRECT_H

#include <mutation_accumulation/statistics/statistics/moment_statistics_population.h> // monte_carlo::Moment_Statistics_Population
#include <mutation_accumulation/utility/labels.h> // labels::create_label

/*************************************************************************/

namespace monte_carlo {

    /** 
     * print moment of population probability distribution versus time for all sub-populations and all species 
     */
    template<class Moment_type, class Configuration_type>
    const double print_moment_populations(
    const Moment_Statistics_Population<Moment_type, Configuration_type> &statistics) {

        typedef Moment_Statistics_Population<Moment_type, Configuration_type> statistics_type;
        typedef typename statistics_type::Results_t Results_t;

        const array::Array3D<Results_t> array3D_results = statistics.get_results_so_far();

        for (int pop = 0; pop < array3D_results.get_dim0(); pop++) {

            const std::string label = labels::create_label(Moment_type());
            const std::string filename = label + "PopulationVersusTime__pop" + boost::lexical_cast<std::string > (pop) + ".dat";
            boost::shared_ptr<std::ofstream> ofstream_ptr = monte_carlo::open_file_for_output(filename);

            for (int node = 0; node < array3D_results.get_dim2(); node++) {

                typedef typename statistics_type::Configuration_t::time_t time_type;
                const time_type tt = statistics.time_grid().at(Node(node));
                *ofstream_ptr << std::setw(30) << std::setprecision(20) << tt;

                for (int spe = 0; spe < array3D_results.get_dim1(); spe++) {

                    const Results_t moment = array3D_results.at(pop, spe, node);
                    *ofstream_ptr << std::setw(30) << std::setprecision(20) << moment;

                }

                *ofstream_ptr << std::endl;

            }
        }
    }
}

#endif	/* MOMENT_POPULATION_DIRECT_H */

