#ifndef BRANCHING_DISCRETE_QUADRATIC_H
#define	BRANCHING_DISCRETE_QUADRATIC_H

#include <boost/assign/list_of.hpp> // boost::assign::list_of()

#include <mutation_accumulation/parameters/parameters_fwd.h> // Symmetry, etc
#include <mutation_accumulation/array/array2D.h> // array::Array2D
#include <mutation_accumulation/random/multinomial_distribution.h> // mutation_accumulation::random::multinomial_distribution
#include <mutation_accumulation/utility/configuration_traits.h> // configuration_types::branching_type etc

#include "configuration.h" // monte_carlo::Configuration
#include "mutation_rates.h" // monte_carlo::MutationRates

/*************************************************************************/

namespace monte_carlo {

    namespace branching_discrete_quadratic {

        /* time is discrete */
        typedef int time_type;

        /**
         * contains state of a branching process in discrete time
         * includes the possibility that mutations occur simultaneously in both daughter stem cells
         */
        template <class population_type>
        class Branching_Discrete_Quadratic : public Configuration<time_type, population_type> {
        public:

            typedef configuration_categories::branching_category category;

        private:

            typedef Configuration<time_type, population_type> base_type;

        private:

            /* probabilities required to implement update_populations_and_time(...) method */
            const MutationRates _uu;
            const Symmetry _ss;
            const SymmetricRenewal _rr;

        private:

            /**
             * calculate probabilities of a cell division lying in reaction category 0 - 5
             */
            const std::vector<double> calculate_categorical_probabilities(const int &spe) const {

                if ((spe >= 0) && (spe < (this->number_species() - 1))) {

                    const double us = _uu.at(spe) * _uu.at(spe);

                    const double p0 = _rr.value() * _ss.value() * (1.0 - 2.0 * _uu.at(spe) + us);
                    const double p1 = _rr.value() * _ss.value() * 2.0 * _uu.at(spe) * (1 - _uu.at(spe));
                    const double p2 = _rr.value() * _ss.value() * us;
                    const double p3 = (1.0 - _ss.value()) * (1.0 - _uu.at(spe));
                    const double p4 = (1.0 - _ss.value()) * _uu.at(spe);
                    const double p5 = (1.0 - _rr.value()) * _ss.value();

                    return boost::assign::list_of(p0) (p1) (p2) (p3) (p4) (p5);

                } else if (spe == (this->number_species() - 1)) {

                    const double p0 = 0.0;
                    const double p1 = 0.0;
                    const double p2 = 0.0;
                    const double p3 = 1.0;
                    const double p4 = 0.0;
                    const double p5 = 0.0;

                    return boost::assign::list_of(p0) (p1) (p2) (p3) (p4) (p5);

                } else {

                    assert(false);
                }
            }

            /**
             * create a matrix to store the random number of times each reaction category occurs
             */
            const array::Array2D<population_type> create_random_matrix(const int &pop, base_generator_type &base_rand_gen) const {

                /* create a matrix to store the random number of times each reaction category occurs for each species */
                const int number_rxn_categories = 6;
                array::Array2D<population_type> R_matrix(this->number_species(), number_rxn_categories);

                /* choose a species */
                for (int spe = 0; spe < this->number_species(); spe++) {

                    /* size of population */
                    const population_type NN = this->get_population(Pop(pop), Spe(spe));

                    /* calculate probabilities of a cell division lying in reaction category 0 - 5 */
                    const std::vector<double> probabilities = calculate_categorical_probabilities(spe);

                    /* create multinomial_distribution object */
                    mutation_accumulation::random::multinomial_distribution<population_type, double> mn_rnd(NN, probabilities);

                    /* randomly choose how many times each reaction category occurs in this sub-population */
                    std::vector<population_type> random_vector = mn_rnd(base_rand_gen);

                    /* error catching */
                    assert(random_vector.size() == R_matrix.get_dim1());

                    /* assign random vector to corresponding row of random matrix \n
                     * note: non-trivial to implement this as: \n
                     * R_matrix.slice0(spe) = random_vector; */
                    for (int cat = 0; cat < R_matrix.get_dim1(); cat++)
                        R_matrix.at(spe, cat) = random_vector.at(cat);

                }

                return R_matrix;
            }

            /**
             * calculate contribution to spe sub-population from divisions of type-(spe-1) cells
             */
            const population_type contribution_from_prior_species(const int &spe, const array::Array2D<population_type> &R_matrix) const {

                assert(spe > 0);

                return (R_matrix.at(spe - 1, 1) + 2.0*R_matrix.at(spe - 1, 2) + R_matrix.at(spe - 1, 4));

            }

            /**
             * calculate contribution to spe sub-population from divisions of type-spe cells
             */
            const population_type contribution_from_current_species(const int &spe, const array::Array2D<population_type> &R_matrix) const {

                return (-R_matrix.at(spe, 0) + R_matrix.at(spe, 2) + R_matrix.at(spe, 4) + R_matrix.at(spe, 5));

            }

            /**
             * update a single independent sub-population
             */
            void update_sub_population(const int &pop, base_generator_type &base_rand_gen) {

                /* create a matrix to store the random number of times each reaction category occurs for each species */
                array::Array2D<population_type> R_matrix_pop = create_random_matrix(pop, base_rand_gen);

                /* update population sizes */
                {
                    /* species = 0 */
                    {
                        const int spe = 0;
                        const population_type old_population_size = this->get_population(Pop(pop), Spe(spe));
                        const population_type new_population_size = old_population_size - contribution_from_current_species(spe, R_matrix_pop);
                        set_population(Pop(pop), Spe(spe), new_population_size);
                    }

                    /* species = 1 ... (number_species - 1) */
                    for (int spe = 1; spe < this->number_species(); spe++) {
                        const population_type old_population_size = this->get_population(Pop(pop), Spe(spe));
                        const population_type new_population_size = old_population_size + contribution_from_prior_species(spe, R_matrix_pop) - contribution_from_current_species(spe, R_matrix_pop);
                        set_population(Pop(pop), Spe(spe), new_population_size);
                    }

                }

            }

            /**
             * update all sub-populations independently 
             */
            void update_populations(base_generator_type &base_rand_gen) {

                /* update each sub-population independently */
                for (int pop = 0; pop < this->number_sub_pops(); pop++)
                    update_sub_population(pop, base_rand_gen);

            }

            /**
             * update time
             */
            void update_time() {

                set_time(this->get_time() + static_cast<time_type> (1));

            }

            /**
             * update populations and time
             */
            virtual void update_populations_and_time(base_generator_type &base_rand_gen) {

                /* update all sub populations */
                update_populations(base_rand_gen);

                /* update time */
                update_time();

            }

        public:

            /**
             * constructor
             */
            explicit Branching_Discrete_Quadratic(
                    const Population2D<population_type> &population2D,
                    const MutationRates &uu,
                    const Symmetry &ss,
                    const SymmetricRenewal &rr,
                    const Uniform_Time_Grid<time_type> &time_grid)
            : base_type(population2D, time_grid), _uu(uu), _ss(ss), _rr(rr) {

                assert(uu.size() == (population2D.number_species() - 1));

            }

        };


    }

    using branching_discrete_quadratic::Branching_Discrete_Quadratic;

}

#endif	/* BRANCHING_DISCRETE_QUADRATIC_H */

