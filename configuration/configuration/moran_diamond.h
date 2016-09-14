#ifndef MORAN_DIAMOND_H
#define	MORAN_DIAMOND_H

#include <boost/random/uniform_real_distribution.hpp> // uniform_real_distribution
#include <boost/random/variate_generator.hpp> // variate_generator
#include <boost/random/exponential_distribution.hpp> // exponential_distribution

#include <mutation_accumulation/parameters/parameters_fwd.h> // Symmetry, etc
#include <mutation_accumulation/array/array4D.h> // array::Array4D
#include <mutation_accumulation/utility/configuration_traits.h> // configuration_types::branching_type etc

//#define DEBUG_MORAN_DIAMOND
#ifdef DEBUG_MORAN_DIAMOND
#include <mutation_accumulation/array/array4D_convenience_functions.h> // array::store_4D_array
#endif

#include "configuration.h" // monte_carlo::Configuration
#include "mutation_rates.h" // monte_carlo::MutationRates

/*************************************************************************/

namespace monte_carlo {

    namespace moran_diamond {

        /* time is cts in a moran model */
        typedef double time_type;

        // indices for each species
        const int spe_0 = 0; // AB
        const int spe_a = 1; // aB
        const int spe_b = 2; // Ab
        const int spe_ab = 3; // AB

        // indices for each mutation rate
        const int mut_0a = 0; // mutation rate from AB to aB
        const int mut_0b = 1; // mutation rate from AB to Ab
        const int mut_a = 2; // mutation rate from aB to ab
        const int mut_b = 3; // mutation rate from Ab to ab

        /**
         * two pathways to a double mutant stem cell \n
         * assumes that Population2D obeys the convention: N_AB N_aB N_Ab N_ab \n
         * assumes that MutationRates obeys the convention: u_0a u_0b u_a u_b \n
         * moran process in cts time \n
         */
        template <class population_type>
        class Moran_Diamond : public Configuration<time_type, population_type> {
        public:

            typedef configuration_categories::moran_category category;

        private:

            typedef Configuration<time_type, population_type> base_type;

        private:

            /* probabilities required to implement update_populations_and_time(...) method */
            const MutationRates _uu;
            const Symmetry _ss;

            const long long int _NN; // total population size (constant in time)

        private:

            /**
             * fetch current population size of stage-spe cells in sub-population pop
             */
            const population_type nn(const int &pop, const int &spe) const {

                return this->get_population(Pop(pop), Spe(spe));
            }

            /**
             * calculate "marginal" transition rates given the current state
             */
            const array::Array4D<double> compute_transition_rates() const {

                // all elements initialized to zero
                array::Array4D<double> ww(
                        this->number_sub_pops(),
                        3,
                        this->number_sub_pops(),
                        4);

                const double ss = _ss.value();

                /* 
                 * jk -> il means that a stage-k cell in the jth sub-population
                 * is converted to a stage-l cell in the ith sub-population
                 */
                for (int pop_dec = 0; pop_dec < this->number_sub_pops(); pop_dec++)
                    for (int pop_inc = 0; pop_inc < this->number_sub_pops(); pop_inc++) {

                        { // j0 -> ix
                            { // j0 -> i0
                                if (pop_dec != pop_inc) {
                                    ww.at(pop_dec, spe_0, pop_inc, spe_0) = 0.5 * ss * nn(pop_inc, spe_0) * (1.0 - 2.0 * (_uu.at(mut_0a) + _uu.at(mut_0b))) * nn(pop_dec, spe_0) / (double) _NN;
                                }
                            }

                            { // j0 -> ia
                                const double lambda_sm = ss * nn(pop_inc, spe_0) * _uu.at(mut_0a) * nn(pop_dec, spe_0) / (double) _NN;
                                double lambda_am = 0.0;
                                if (pop_dec == pop_inc)
                                    lambda_am = (1.0 - ss) * nn(pop_inc, spe_0) * _uu.at(mut_0a);
                                const double lambda_s = 0.5 * ss * nn(pop_inc, spe_a) * (1.0 - 2.0 * _uu.at(mut_a)) * nn(pop_dec, spe_0) / (double) _NN;
                                ww.at(pop_dec, spe_0, pop_inc, spe_a) = lambda_sm + lambda_am + lambda_s;
                            }

                            { // j0 -> ib
                                const double lambda_sm = ss * nn(pop_inc, spe_0) * _uu.at(mut_0b) * nn(pop_dec, spe_0) / (double) _NN;
                                double lambda_am = 0.0;
                                if (pop_dec == pop_inc)
                                    lambda_am = (1.0 - ss) * nn(pop_inc, spe_0) * _uu.at(mut_0b);
                                const double lambda_s = 0.5 * ss * nn(pop_inc, spe_b) * (1.0 - 2.0 * _uu.at(mut_b)) * nn(pop_dec, spe_0) / (double) _NN;
                                ww.at(pop_dec, spe_0, pop_inc, spe_b) = lambda_sm + lambda_am + lambda_s;
                            }

                            { // j0 -> i,ab
                                const double lambda_sm_a = ss * nn(pop_inc, spe_a) * _uu.at(mut_a) * nn(pop_dec, spe_0) / (double) _NN;
                                const double lambda_sm_b = ss * nn(pop_inc, spe_b) * _uu.at(mut_b) * nn(pop_dec, spe_0) / (double) _NN;
                                ww.at(pop_dec, spe_0, pop_inc, spe_ab) = lambda_sm_a + lambda_sm_b;
                            }
                        }
                        { // ja -> ix
                            { // ja -> i0
                                ww.at(pop_dec, spe_a, pop_inc, spe_0) = 0.5 * ss * nn(pop_inc, spe_0) * (1.0 - 2.0 * (_uu.at(mut_0a) + _uu.at(mut_0b))) * nn(pop_dec, spe_a) / (double) _NN;
                            }

                            { // ja -> ia
                                if (pop_dec != pop_inc) {
                                    const double lambda_sm = ss * nn(pop_inc, spe_0) * _uu.at(mut_0a) * nn(pop_dec, spe_a) / (double) _NN;
                                    const double lambda_s = 0.5 * ss * nn(pop_inc, spe_a) * (1.0 - 2.0 * _uu.at(mut_a)) * nn(pop_dec, spe_a) / (double) _NN;
                                    ww.at(pop_dec, spe_a, pop_inc, spe_a) = lambda_sm + lambda_s;
                                }
                            }

                            { // ja -> ib
                                const double lambda_sm = ss * nn(pop_inc, spe_0) * _uu.at(mut_0b) * nn(pop_dec, spe_a) / (double) _NN;
                                const double lambda_s = 0.5 * ss * nn(pop_inc, spe_b) * (1.0 - 2.0 * _uu.at(mut_b)) * nn(pop_dec, spe_a) / (double) _NN;
                                ww.at(pop_dec, spe_a, pop_inc, spe_b) = lambda_sm + lambda_s;
                            }

                            { // ja -> i,ab
                                const double lambda_sm_a = ss * nn(pop_inc, spe_a) * _uu.at(mut_a) * nn(pop_dec, spe_a) / (double) _NN;
                                const double lambda_sm_b = ss * nn(pop_inc, spe_b) * _uu.at(mut_b) * nn(pop_dec, spe_a) / (double) _NN;
                                double lambda_am = 0.0;
                                if (pop_dec == pop_inc)
                                    lambda_am = (1.0 - ss) * nn(pop_inc, spe_a) * _uu.at(mut_a);
                                ww.at(pop_dec, spe_a, pop_inc, spe_ab) = lambda_sm_a + lambda_am + lambda_sm_b;
                            }
                        }
                        { // jb -> ix
                            { // jb -> i0
                                ww.at(pop_dec, spe_b, pop_inc, spe_0) = 0.5 * ss * nn(pop_inc, spe_0) * (1.0 - 2.0 * (_uu.at(mut_0a) + _uu.at(mut_0b))) * nn(pop_dec, spe_b) / (double) _NN;
                            }

                            { // jb -> ia
                                const double lambda_sm = ss * nn(pop_inc, spe_0) * _uu.at(mut_0a) * nn(pop_dec, spe_b) / (double) _NN;
                                const double lambda_s = 0.5 * ss * nn(pop_inc, spe_a) * (1.0 - 2.0 * _uu.at(mut_a)) * nn(pop_dec, spe_b) / (double) _NN;
                                ww.at(pop_dec, spe_b, pop_inc, spe_a) = lambda_sm + lambda_s;
                            }

                            { // jb -> ib
                                if (pop_dec != pop_inc) {
                                    const double lambda_sm = ss * nn(pop_inc, spe_0) * _uu.at(mut_0b) * nn(pop_dec, spe_b) / (double) _NN;
                                    const double lambda_s = 0.5 * ss * nn(pop_inc, spe_b) * (1.0 - 2.0 * _uu.at(mut_b)) * nn(pop_dec, spe_b) / (double) _NN;
                                    ww.at(pop_dec, spe_b, pop_inc, spe_b) = lambda_sm + lambda_s;
                                }
                            }

                            { // jb -> i,ab
                                const double lambda_sm_a = ss * nn(pop_inc, spe_a) * _uu.at(mut_a) * nn(pop_dec, spe_b) / (double) _NN;
                                const double lambda_sm_b = ss * nn(pop_inc, spe_b) * _uu.at(mut_b) * nn(pop_dec, spe_b) / (double) _NN;
                                double lambda_am = 0.0;
                                if (pop_dec == pop_inc)
                                    lambda_am = (1.0 - ss) * nn(pop_inc, spe_b) * _uu.at(mut_b);
                                ww.at(pop_dec, spe_b, pop_inc, spe_ab) = lambda_sm_a + lambda_am + lambda_sm_b;
                            }
                        }
                    }


                return ww;

            }

            /**
             * change number of spe-type mutants in sub-population pop by -1 or +1
             */
            void update_populations(const int &pop, const int &spe, const int &update_value) {

                assert((update_value == -1) || (update_value == +1));

                // explicitly make the lookup of following base method template-parameter-dependent by calling it through this->
                const population_type old_population_size = this->get_population(Pop(pop), Spe(spe));
                const population_type new_population_size = old_population_size + static_cast<population_type> (update_value);
                // this-> not needed in the following base method call because it has arguments that depend on a template parameter
                // http://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
                // http://eli.thegreenplace.net/2012/02/06/dependent-name-lookup-for-c-templates/
                set_population(Pop(pop), Spe(spe), new_population_size);

            }

            /**
             * update populations
             */
            void update_populations_gillespie(
                    base_generator_type &base_rand_gen,
                    const array::Array4D<double> &ww,
                    const double &ww_total) {

                /* randomly choose next reaction using inverse CDF method */
                typedef boost::random::uniform_real_distribution<> uniform_generator_type;
                typedef boost::random::variate_generator<base_generator_type&, uniform_generator_type> uniform_variate_type;
                uniform_variate_type uniform_random_number(base_rand_gen, uniform_generator_type(0.0, 1.0));
                const double random_fraction_ww_total = uniform_random_number() * ww_total;
                int pop_dec, spe_dec, pop_inc, spe_inc;
                if (!ww.cumulative_sum(random_fraction_ww_total, pop_dec, spe_dec, pop_inc, spe_inc)) {
                    std::cerr << "could not choose next reaction" << std::endl;
                    std::cerr << "random_fraction_ww_total = " << random_fraction_ww_total << std::endl;
                    std::cerr << "ww_total = " << ww_total << std::endl;
                    if (ww_total < 1e-10)
                        std::cerr << "process could be stuck in absorbing state" << std::endl;
                    print_debug_info(*this);
                    assert(false);
                }

                /* execute chosen reaction */
                update_populations(pop_dec, spe_dec, -1);
                update_populations(pop_inc, spe_inc, +1);

            }

            /**
             * update time
             */
            void update_time_gillespie(base_generator_type &base_rand_gen, const double &ww_total) {

                /* update time using an exponentially distributed inter-event time */
                typedef boost::random::exponential_distribution<> exponential_generator_type;
                typedef boost::random::variate_generator<base_generator_type&, exponential_generator_type> exponential_variate_type;
                exponential_variate_type exponential_random_number(base_rand_gen, exponential_generator_type(ww_total));

                // this-> makes argument to set_time template-parameter-dependent 
                set_time(this->get_time() + static_cast<time_type> (exponential_random_number()));

            }

            /**
             * update populations and time
             */
            virtual void update_populations_and_time(base_generator_type &base_rand_gen) {

                /* calculate "marginal" transition rates
                 * ie rates at which system changes state 
                 */
                const array::Array4D<double> ww = compute_transition_rates();

#ifdef DEBUG_MORAN_DIAMOND
                /* inspect ww array */
                array::store_4D_array(ww, "ww.dat");
                std::cout << "ww = " << ww.sum() << std::endl;
#endif

                /* calculate rate at which next event occurs */
                const double ww_total = ww.sum();

                /* update populations */
                update_populations_gillespie(base_rand_gen, ww, ww_total);

                /* update time */
                update_time_gillespie(base_rand_gen, ww_total);

            }

        public:

            /**
             * constructor
             */
            explicit Moran_Diamond(
                    const Population2D<population_type> &population2D,
                    const MutationRates &uu,
                    const Symmetry &ss,
                    const Uniform_Time_Grid<time_type> &time_grid)
            : base_type(population2D, time_grid), _uu(uu), _ss(ss), _NN(population2D.sum()) {

                assert(population2D.number_species() == 4);
                assert(uu.size() == 4);
            }

        };


    }

    using moran_diamond::Moran_Diamond;

}

#endif	/* MORAN_DIAMOND_H */

