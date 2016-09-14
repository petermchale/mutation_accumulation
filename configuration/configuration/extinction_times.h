#ifndef EXTINCTION_TIMES_H
#define	EXTINCTION_TIMES_H

#include <vector> // std::vector

#include "population2D.h" // monte_carlo::Population2D

/*************************************************************************/

namespace monte_carlo {


    /**
     * store times at which each sub-population extinguished
     */
    template <class time_type, class population_type>
    class ExtinctionTimes {

    private:

        /** 
         * _extinction_times.at(pop) is time at which sub-population pop extinguished \n
         * value of -1 indicates that extinction has not yet occurred \n
         **/
        std::vector<time_type> _extinction_times;


    public:
        
        /**
         * custom constructor
         */
        explicit ExtinctionTimes(const int &number_sub_pops) {

            _extinction_times = std::vector<time_type>(number_sub_pops, static_cast<time_type>(-1));

        }

        
        /**
         * update extinction times
         */
        void update(const time_type &tt, const Population2D<population_type> &population2D) {

            for (int pop = 0; pop < population2D.number_sub_pops(); pop++) // choose a sub-population ...

                if (_extinction_times.at(pop) < static_cast<time_type>(0)) // if an extinction has not been recorded ...

                    if (sub_pop_size(population2D, Pop(pop)) == static_cast<population_type>(0)) // if an extinction has actually occurred ...

                        _extinction_times.at(pop) = tt; // record the extinction

        }


        /**
         * get current extinction times
         */
        const time_type operator()(const Pop &pop) const {

            return _extinction_times.at(pop.value());
        }
        

    }; 

    
} 

#endif	/* EXTINCTION_TIMES_H */

