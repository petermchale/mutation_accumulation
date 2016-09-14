#ifndef MORAN3_H
#define	MORAN3_H

#include <mutation_accumulation/utility/configuration_traits.h> // configuration_types::branching_type etc

#include "population2D.h" // monte_carlo::Population2D
#include "mutation_rates.h" // monte_carlo::MutationRates
#include "configuration_interface.h" // monte_carlo::Configuration_Interface
#include "moran4.h" // monte_carlo::Moran4 

//#define DEBUG_MORAN3_H

/*************************************************************************/

namespace monte_carlo {

    namespace moran3 {

        namespace moran3_utility {

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
         * accumulation of two mutations (3 stages) in a Moran process \n
         * requires three species and two mutation rates\n
         * delegates calculations to Moran4 (see "delegation pattern" in wikipedia and p.352 of Fowler)
         */
        template <class population_type>
        class Moran3 : public Configuration_Interface<time_type, population_type> {
        public:

            typedef configuration_categories::moran_category category;

        private:

            typedef Population2D<population_type> population2D_type;

        private:

            Moran4<population_type> _moran4;

        private:

            /**
             * create a 4-species Population2D object from given 3-species Population2D object\n
             * static function removes danger of accidentally referring to the nascent Moran3 object's as-yet-uninitialized data members (item 9)
             */
            static const population2D_type convert_3Species_to_4Species(const population2D_type &population2D_3species) {

                assert(population2D_3species.number_species() == 3);

                population2D_type population2D_4species(population2D_3species.number_sub_pops(), 4);

                for (int pop = 0; pop < population2D_3species.number_sub_pops(); pop++) {
                    for (int spe = 0; spe < population2D_3species.number_species(); spe++)
                        population2D_4species.at(pop, spe) = population2D_3species.at(pop, spe);
                    population2D_4species.at(pop, population2D_3species.number_species()) = static_cast<population_type> (0);
                }

                return population2D_4species;
            }

            /**
             * create a 4-species MutationRates object from a 3-species MutationRates object\n
             * static function removes danger of accidentally referring to the nascent Moran3 object's as-yet-uninitialized data members (item 9)
             */
            static const MutationRates convert_3Species_to_4Species(const MutationRates &mutationRates_3species) {

                assert(mutationRates_3species.size() == 2);

                MutationRates mutationRates_4species(3);

                for (int spe = 0; spe < mutationRates_3species.size(); spe++)
                    mutationRates_4species.at(spe) = mutationRates_3species.at(spe);
                mutationRates_4species.at(mutationRates_3species.size()) = static_cast<MutationRates::data_t> (0);

                return mutationRates_4species;
            }

        public:

            /**
             * constructor
             */
            explicit Moran3(
                    const population2D_type &population2D,
                    const MutationRates &uu,
                    const Symmetry &ss,
                    const Uniform_Time_Grid<time_type> &time_grid)
            : _moran4(convert_3Species_to_4Species(population2D), convert_3Species_to_4Species(uu), ss, time_grid) {

#ifdef DEBUG_MORAN3_H
                moran3_utility::debug_print(convert_3Species_to_4Species(population2D));
                moran3_utility::debug_print(convert_3Species_to_4Species(uu));
#endif

            }

            /**
             * get current time
             */
            virtual const time_type get_time() const {

                return _moran4.get_time();
            }

            /**
             * get number of type-spe cells in sub-population pop
             */
            virtual const population_type get_population(const Pop &pop, const Spe &spe) const {

                return _moran4.get_population(pop, spe);
            }

            /**
             * number of sub-populations
             */
            virtual const int number_sub_pops() const {

                return _moran4.number_sub_pops();
            }

            /**
             * number of species
             */
            virtual const int number_species() const {

                return 3;
            }

            /**
             * number of filled nodes
             */
            virtual const int number_filled_nodes() const {

                return _moran4.number_filled_nodes();
            }

            /**
             * get time at which sub-population pop extinguished
             */
            virtual const time_type extinction_times(const Pop & pop) const {

                return _moran4.extinction_times(pop);
            }

            /**
             * get time at which first spe-type cell appears in sub-population pop
             */
            virtual const time_type mutation_times(const Pop &pop, const Spe & spe) const {

                return _moran4.mutation_times(pop, spe);
            }

            /**
             * get path time \n
             */
            virtual const time_type get_path_time(const Node &node) const {

                return _moran4.get_path_time(node);
            }

            /**
             * get time of last node in path\n
             */
            virtual const time_type get_last_node_time() const {

                return _moran4.get_last_node_time();
            }

            /**
             * get path population \n
             */
            virtual const population_type get_path_population(const Pop &pop, const Spe &spe, const Node &node) const {

                return _moran4.get_path_population(pop, spe, node);
            }

            /**
             * transition to next state of stochastic process
             */
            virtual void transition(base_generator_type &base_rand_gen) {

                _moran4.transition(base_rand_gen);
            }

            /**
             * determine if end of trajectory has been reached
             */
            virtual const bool end() const {

                return _moran4.end();

            }


        };


    }

    using moran3::Moran3;

}

#endif	/* MORAN3_H */

