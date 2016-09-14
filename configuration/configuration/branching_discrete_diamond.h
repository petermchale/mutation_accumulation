#ifndef BRANCHING_DISCRETE_DIAMOND_H
#define	BRANCHING_DISCRETE_DIAMOND_H

#include <boost/assign/list_of.hpp> // boost::assign::list_of()

#include <mutation_accumulation/random/multinomial_distribution.h> // mutation_accumulation::random::multinomial_distribution
#include <mutation_accumulation/utility/configuration_traits.h> // configuration_categories::branching_category

#include "configuration.h" // monte_carlo::Configuration
#include "route_along_diamond.h" // Diamond_Routes
#include "mutation_rates.h" // monte_carlo::MutationRates

/*************************************************************************/

namespace monte_carlo {

    namespace branching_discrete_diamond {

        /* time is discrete */
        typedef int time_type;

        // indices for each species
        const int spe_0 = 0; // AB
        const int spe_a = 1; // aB
        const int spe_b = 2; // Ab
        const int spe_ab = 3; // AB

        // indices for multiple species 
        const std::vector<int> wt_singleMutant_indices = boost::assign::list_of(spe_0)(spe_a) (spe_b);
        const std::vector<int> singleMutant_indices = boost::assign::list_of(spe_a) (spe_b);

        // indices for each mutation rate
        const int mut_0a = 0; // mutation rate from AB to aB
        const int mut_0b = 1; // mutation rate from AB to Ab
        const int mut_a = 2; // mutation rate from aB to ab
        const int mut_b = 3; // mutation rate from Ab to ab

        /**
         * two pathways to a double mutant stem cell \n
         * assumes that Population2D obeys the convention: N_AB N_aB N_Ab N_ab \n
         * assumes that MutationRates obeys the convention: u_0a u_0b u_a u_b \n
         * branching process in discrete time \n
         */
        template <class population_type>
        class Branching_Discrete_Diamond : public Configuration<time_type, population_type> {
        public:

            typedef configuration_categories::branching_category category;

        private:

            typedef Configuration<time_type, population_type> base_type;

        private:

            /* probabilities required to implement update_populations_and_time(...) method */
            const MutationRates _uu;
            const Symmetry _ss;
            const SymmetricRenewal _rr;

            Diamond_Routes<population_type> _diamond_routes;

        private:

            /**
             * determine if spe is wild-type species
             */
            const bool wild_type(const int &spe) const {

                return (spe == spe_0);

            }

            /**
             * determine if spe is a single-mutant species
             */
            const bool single_mutant(const int &spe) const {

                return ((spe == spe_a) || (spe == spe_b));

            }

            /**
             * calculate mutation rate from single-mutant to double-mutant 
             */
            const double calculate_uu_single_mutant(const int &spe) const {

                if (spe == spe_a)
                    return _uu.at(mut_a);
                else if (spe == spe_b)
                    return _uu.at(mut_b);
                else
                    assert(false);

            }

            /**
             * calculate probabilities of a cell division of type-spe lying in reaction category-cat
             */
            const std::vector<double> calculate_categorical_probabilities(const int &spe) const {

                if (wild_type(spe)) {

                    // unpack relevant mutation rates
                    const double uu_0a = _uu.at(mut_0a);
                    const double uu_0b = _uu.at(mut_0b);

                    const double p0 = _rr.value() * _ss.value() * (1.0 - 2.0 * (uu_0a + uu_0b));
                    const double p1 = _rr.value() * _ss.value() * 2.0 * uu_0a;
                    const double p2 = _rr.value() * _ss.value() * 2.0 * uu_0b;
                    const double p3 = (1.0 - _ss.value()) * (1.0 - (uu_0a + uu_0b));
                    const double p4 = (1.0 - _ss.value()) * uu_0a;
                    const double p5 = (1.0 - _ss.value()) * uu_0b;
                    const double p6 = (1.0 - _rr.value()) * _ss.value();

                    return boost::assign::list_of(p0) (p1) (p2) (p3) (p4) (p5) (p6);

                } else if (single_mutant(spe)) {

                    // unpack relevant mutation rate
                    const double uu_single_mutant = calculate_uu_single_mutant(spe);

                    const double p0 = _rr.value() * _ss.value() * (1.0 - 2.0 * uu_single_mutant);
                    const double p1 = _rr.value() * _ss.value() * 2.0 * uu_single_mutant;
                    const double p2 = (1.0 - _ss.value()) * (1.0 - uu_single_mutant);
                    const double p3 = (1.0 - _ss.value()) * uu_single_mutant;
                    const double p4 = (1.0 - _rr.value()) * _ss.value();

                    return boost::assign::list_of(p0) (p1) (p2) (p3) (p4);

                } else {

                    assert(false);
                }
            }

            /**
             * create a vector of vectors to store the random number of times each reaction category occurs for each species (excluding double mutant)
             */
            const std::vector<std::vector<population_type> > draw_zz(const int &pop, base_generator_type &base_rand_gen) const {

                std::vector<std::vector<population_type> > zz(wt_singleMutant_indices.size());

                /* choose a species */
                for (int spe = 0; spe < zz.size(); spe++) {

                    /* size of population */
                    const population_type NN = this->get_population(Pop(pop), Spe(spe));

                    /* calculate probabilities of a cell division lying in each possible reaction category */
                    const std::vector<double> probabilities = calculate_categorical_probabilities(spe);

                    /* create multinomial_distribution object */
                    mutation_accumulation::random::multinomial_distribution<population_type, double> mn_rnd(NN, probabilities);

                    /* randomly choose how many times each reaction category occurs */
                    zz.at(spe) = mn_rnd(base_rand_gen);

                }

                return zz;
            }

            /**
             * mutation of wild-type to single mutant (spe)
             */
            const population_type A0(const int &spe, const std::vector<std::vector<population_type> > &zz) const {

                if (spe == spe_a) {

                    return (zz.at(spe_0).at(1) + zz.at(spe_0).at(4));

                } else if (spe == spe_b) {

                    return (zz.at(spe_0).at(2) + zz.at(spe_0).at(5));

                } else {

                    assert(false);
                }
            }

            /**
             * mutation of single mutant (spe) to double mutant 
             */
            const population_type A1(const int &spe, const std::vector<std::vector<population_type> > &zz) const {

                if (single_mutant(spe)) {

                    return (zz.at(spe).at(1) + zz.at(spe).at(3));

                } else {

                    assert(false);
                }
            }

            /**
             * calculate contribution to spe sub-population from divisions of cells within that sub-population
             */
            const population_type contribution_from_current_species(const int &spe, const std::vector<std::vector<population_type> > &zz) const {

                if (wild_type(spe)) {

                    return (-zz.at(spe).at(0) + zz.at(spe).at(4) + zz.at(spe).at(5) + zz.at(spe).at(6));

                } else if (single_mutant(spe)) {

                    return (-zz.at(spe).at(0) + zz.at(spe).at(3) + zz.at(spe).at(4));

                } else {

                    assert(false);
                }

            }

            /**
             * update a single independent sub-population
             */
            void update_sub_population(const int &pop, base_generator_type &base_rand_gen) {

                /* draw random number of times each reaction category occurs for each species */
                std::vector<std::vector<population_type> > zz = draw_zz(pop, base_rand_gen);

                /* update population sizes */
                {
                    /* wt species */
                    {
                        const population_type old_population_size = this->get_population(Pop(pop), Spe(spe_0));
                        const population_type new_population_size = old_population_size - contribution_from_current_species(spe_0, zz);
                        set_population(Pop(pop), Spe(spe_0), new_population_size);
                    }

                    /* single-mutants */
                    for (int ii = 0; ii < singleMutant_indices.size(); ii++) {
                        const int spe = singleMutant_indices.at(ii);
                        const population_type old_population_size = this->get_population(Pop(pop), Spe(spe));
                        const population_type new_population_size = old_population_size + A0(spe, zz) - contribution_from_current_species(spe, zz);
                        set_population(Pop(pop), Spe(spe), new_population_size);
                    }

                    /* double-mutant */
                    {
                        const population_type old_population_size = this->get_population(Pop(pop), Spe(spe_ab));
                        const population_type new_population_size = old_population_size + A1(spe_a, zz) + A1(spe_b, zz);
                        set_population(Pop(pop), Spe(spe_ab), new_population_size);

                    }

                }

                /* query whether a and/or b generated ab */
                _diamond_routes.update(pop, A1(spe_a, zz), A1(spe_b, zz));
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
            explicit Branching_Discrete_Diamond(
                    const Population2D<population_type> &population2D,
                    const MutationRates &uu,
                    const Symmetry &ss,
                    const SymmetricRenewal &rr,
                    const Uniform_Time_Grid<time_type> &time_grid)
            : base_type(population2D, time_grid), _uu(uu), _ss(ss), _rr(rr), _diamond_routes(population2D.number_sub_pops()) {

                assert(population2D.number_species() == 4);
                assert(uu.size() == 4);
            }

            /**
             * query route a \n
             * generalize to a virtual function (defined in Configuration_Interface) that returns a bit string of length K indicating whether or not each K-1-fold mutant generated first K-fold mutant \n
             */
            const bool a_yielded_ab(const Pop &pop) const {

                return _diamond_routes.a_yielded_ab(pop);

            }

            /**
             * query route b
             */
            const bool b_yielded_ab(const Pop &pop) const {

                return _diamond_routes.b_yielded_ab(pop);

            }


        };


    }

    using branching_discrete_diamond::Branching_Discrete_Diamond;

}

#endif	/* BRANCHING_DISCRETE_DIAMOND_H */

