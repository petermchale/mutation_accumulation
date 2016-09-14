#ifndef MOMENT_STATISTICS_MUTATION_TIME_H
#define	MOMENT_STATISTICS_MUTATION_TIME_H

#include "moment_statistics.h"

/*************************************************************************/

namespace monte_carlo {

    /** 
     * calculate indicated moment of first passage time to Spe-mutant in Pop sub-population 
     */
    template<class Moment_type, class Configuration_type>
    class Moment_Statistics_MutationTime : public Moment_Statistics<Moment_type, Configuration_type> {
    private:

        typedef typename Configuration_type::time_t time_type;

        typedef Moment_Statistics<Moment_type, Configuration_type> base_type;

    public:

        /**
         * constructor \n
         */
        explicit Moment_Statistics_MutationTime(
                const Number_Pop &number_pop,
                const Number_Spe &number_spe,
                const Pop &pop_to_observe,
                const Spe &spe_to_observe,
                const long long int &number_trials_max)
        : base_type(
        number_pop,
        number_spe.value(),
        1,
        number_spe.value(),
        1,
        pop_to_observe,
        spe_to_observe.value(),
        0,
        number_trials_max) {

            /* sample type of Moment_type and time type of Configuration_type should agree */
            BOOST_STATIC_ASSERT((boost::is_same<typename Moment_type::sample_t, typename Configuration_type::time_t>::value));

        }

        /** 
         * dump results of a particular trial
         */
        virtual void dump(const Configuration_type &configuration) {

            for (int pop = 0; pop < configuration.number_sub_pops(); pop++)
                for (int spe = 0; spe < configuration.number_species(); spe++) {
                    if (fate_bool(configuration, Pop(pop), Spe(spe))) {
                        const time_type mutation_time = configuration.mutation_times(Pop(pop), Spe(spe));
                        this->update_moments(Pop(pop), spe, 0, mutation_time);
                    }
                }


            for (int spe = 0; spe < configuration.number_species(); spe++) {
                if (fate_bool(configuration, Spe(spe))) {
                    const time_type mutation_time = mutation_time_whole(configuration, Spe(spe));
                    this->update_moments_whole(spe, 0, mutation_time);
                }
            }


        }

    };



}

#endif	/* MOMENT_STATISTICS_MUTATION_TIME_H */

