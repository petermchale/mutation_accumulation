#ifndef MUTATION_TIMES_H
#define	MUTATION_TIMES_H

#include <mutation_accumulation/array/array2D.h> // array::Array2D

#include "population2D.h" // monte_carlo::Population2D

/*************************************************************************/

namespace monte_carlo {

    /**
     * store times at which first spe-type cell appears in sub-population pop
     */
    template <class time_type, class population_type>
    class MutationTimes {
    private:

        typedef Population2D<population_type> population2D_type;

    private:

        /** 
         * _mutation_times(pop, spe) is time at which first spe-type cell appears in sub-population pop \n
         * value of -1 indicates that mutation has not yet occurred \n
         **/
        array::Array2D<time_type> _mutation_times;

        // boolean variables indicating whether corresponding mutations have occurred or not 
        array::Array2D<bool> _mutations_occurred;


    private:

        /**
         * implement update mutation times
         */
        void implement_update(const time_type &tt, const population2D_type &population2D) {

            for (int pop = 0; pop < population2D.number_sub_pops(); pop++)
                for (int spe = 0; spe < population2D.number_species(); spe++)
                    if (!_mutations_occurred.at(pop, spe)) // if the mutation has not been recorded ...
                        if (population2D.at(pop, spe) > static_cast<population_type> (0)) // if the mutation has actually occurred ...
                        {
                            _mutation_times.at(pop, spe) = tt; // record the mutation time
                            _mutations_occurred.at(pop, spe) = true; // indicate that mutation has occurred
                        }
        }

    public:

        /**
         * custom constructor
         */
        explicit MutationTimes(const population2D_type & population2D) {

            _mutation_times = array::Array2D<time_type > (population2D.number_sub_pops(), population2D.number_species(), static_cast<time_type> (-1));
            _mutations_occurred = array::Array2D<bool > (population2D.number_sub_pops(), population2D.number_species(), false);

            implement_update(static_cast<time_type> (0), population2D);

        }

        /**
         * update mutation times
         */
        void update(const time_type &tt, const population2D_type & population2D) {

            implement_update(tt, population2D);

        }

        /**
         * get current mutation times
         */
        const time_type operator()(const Pop &pop, const Spe & spe) const {

            return _mutation_times.at(pop.value(), spe.value());
        }


    };


}

#endif	/* MUTATION_TIMES_H */

