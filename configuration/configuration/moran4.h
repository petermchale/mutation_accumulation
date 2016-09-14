#ifndef MORAN4_H
#define	MORAN4_H

#include <mutation_accumulation/utility/configuration_traits.h> // configuration_types::branching_type etc

#include "population2D.h" // monte_carlo::Population2D
#include "mutation_rates.h" // monte_carlo::MutationRates
#include "configuration_interface.h" // monte_carlo::Configuration_Interface
#include "moran5.h" // monte_carlo::Moran5

//#define DEBUG_MORAN4

/*************************************************************************/

namespace monte_carlo {

    namespace moran4 {

        namespace moran4_utility {

            template <class population_type>
            const void debug_print(const Population2D<population_type> &population2D) {

                std::cout << "population2D = " << std::endl;
                for (int pop = 0; pop < population2D.number_sub_pops(); pop++) {
                    for (int spe = 0; spe < population2D.number_species(); spe++)
                        std::cout << "<" << population2D.at(pop, spe) << ">" << " ";
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }

            const void debug_print(const MutationRates &mutationRates) {

                std::cout << "mutationRates = " << std::endl;
                for (int spe = 0; spe < mutationRates.size(); spe++)
                    std::cout << "<" << mutationRates.at(spe) << ">" << " ";
                std::cout << std::endl;
            }

        }

        /* time is cts in a moran model */
        typedef double time_type;

        /**
         * accumulation of three mutations (4 stages) in a Moran process \n
         * requires three mutation rates\n
         * delegates calculations to Moran5 (see "delegation pattern" in wikipedia and p.352 of Fowler)
         */
        template <class population_type>
        class Moran4 : public Configuration_Interface<time_type, population_type> {
        public:

            typedef configuration_categories::moran_category category;

        private:

            typedef Population2D<population_type> population2D_type;

        private:

            Moran5<population_type> _moran5;

        private:

            /**
             * create a 5-species Population2D object from given 4-species Population2D object\n
             * static function removes danger of accidentally referring to the nascent Moran4 object's as-yet-uninitialized data members (item 9)
             */
            static const population2D_type convert_4Species_to_5Species(const population2D_type &population2D_4species) {

                assert(population2D_4species.number_species() == 4);

                population2D_type population2D_5species(population2D_4species.number_sub_pops(), 5);

                for (int pop = 0; pop < population2D_4species.number_sub_pops(); pop++) {
                    for (int spe = 0; spe < population2D_4species.number_species(); spe++)
                        population2D_5species.at(pop, spe) = population2D_4species.at(pop, spe);
                    population2D_5species.at(pop, population2D_4species.number_species()) = static_cast<population_type> (0);
                }

                return population2D_5species;
            }

            /**
             * create a 5-species MutationRates object from a 4-species MutationRates object\n
             * static function removes danger of accidentally referring to the nascent Moran4 object's as-yet-uninitialized data members (item 9)
             */
            static const MutationRates convert_4Species_to_5Species(const MutationRates &mutationRates_4species) {

                assert(mutationRates_4species.size() == 3);

                MutationRates mutationRates_5species(4);

                for (int spe = 0; spe < mutationRates_4species.size(); spe++)
                    mutationRates_5species.at(spe) = mutationRates_4species.at(spe);
                mutationRates_5species.at(mutationRates_4species.size()) = static_cast<MutationRates::data_t> (0);

                return mutationRates_5species;
            }

        public:

            /**
             * constructor
             */
            explicit Moran4(
                    const population2D_type &population2D,
                    const MutationRates &uu,
                    const Symmetry &ss,
                    const Uniform_Time_Grid<time_type> &time_grid)
            : _moran5(convert_4Species_to_5Species(population2D), convert_4Species_to_5Species(uu), ss, time_grid) {

#ifdef DEBUG_MORAN4
                moran4_utility::debug_print(convert_4Species_to_5Species(population2D));
                moran4_utility::debug_print(convert_4Species_to_5Species(uu));
#endif

            }

            /**
             * get current time
             */
            virtual const time_type get_time() const {

                return _moran5.get_time();
            }

            /**
             * get number of type-spe cells in sub-population pop
             */
            virtual const population_type get_population(const Pop &pop, const Spe &spe) const {

                return _moran5.get_population(pop, spe);
            }

            /**
             * number of sub-populations
             */
            virtual const int number_sub_pops() const {

                return _moran5.number_sub_pops();
            }

            /**
             * number of species
             */
            virtual const int number_species() const {

                return 4;
            }

            /**
             * number of filled nodes
             */
            virtual const int number_filled_nodes() const {

                return _moran5.number_filled_nodes();
            }

            /**
             * get time at which sub-population pop extinguished
             */
            virtual const time_type extinction_times(const Pop & pop) const {

                return _moran5.extinction_times(pop);
            }

            /**
             * get time at which first spe-type cell appears in sub-population pop
             */
            virtual const time_type mutation_times(const Pop &pop, const Spe & spe) const {

                return _moran5.mutation_times(pop, spe);
            }

            /**
             * get path time \n
             */
            virtual const time_type get_path_time(const Node &node) const {

                return _moran5.get_path_time(node);
            }

            /**
             * get time of last node in path\n
             */
            virtual const time_type get_last_node_time() const {

                return _moran5.get_last_node_time();
            }

            /**
             * get path population \n
             */
            virtual const population_type get_path_population(const Pop &pop, const Spe &spe, const Node &node) const {

                return _moran5.get_path_population(pop, spe, node);
            }

            /**
             * transition to next state of stochastic process
             */
            virtual void transition(base_generator_type &base_rand_gen) {

                _moran5.transition(base_rand_gen);
            }

            /**
             * determine if end of trajectory has been reached
             */
            virtual const bool end() const {

                return _moran5.end();

            }


        };


    }

    using moran4::Moran4;

}

#endif	/* MORAN4_H */

