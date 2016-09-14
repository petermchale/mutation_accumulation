#ifndef PATH_POLICY_H
#define	PATH_POLICY_H

#include <mutation_accumulation/configuration/utilities/fpt_mutation_lastSpecies.h> // monte_carlo::all_last_species_occurred
#include <mutation_accumulation/configuration/utilities/fpt_extinction.h> // monte_carlo::extinguished
#include <mutation_accumulation/configuration/utilities/fate.h> // monte_carlo::eventualFateOccurred_lastSpecies_allSubPops

#include <mutation_accumulation/configuration/configuration/random_fwd.h> // monte_carlo::base_generator_type

/*************************************************************************/

namespace monte_carlo {

    /**
     * abstract base class that evaluates stochastic trajectory\n
     * \n
     * Template Method Pattern via the Non-Virtual Interface Idiom (Item 35)\n
     * \n
     * may use this class polymorphically (virtual destructor) \n
     * may not instantiate (ctor is protected)\n
     */
    template<class Configuration_type>
    class Path_Policy_Base {
    public:

        typedef Configuration_type Configuration_t;

    private:

        typedef typename Configuration_type::time_t time_type; // Item 42
        typedef typename Configuration_type::population_t population_type;

    private:
        
        std::vector<base_generator_type> _successful_PRNG_states;

    private:

        /** 
         * return true if path is successful
         */
        virtual const bool success(const Configuration_type &configuration) const = 0;


    protected:

        /** 
         * protected ctor prevents instantiation of base class 
         * 
         */
        explicit Path_Policy_Base() : _successful_PRNG_states(std::vector<base_generator_type>()) {

            /* check template parameter types to supplement "duck typing" */
            BOOST_STATIC_ASSERT((boost::is_base_of<Configuration_Interface<time_type, population_type>, Configuration_type>::value));

            /* check size of vector is zero */
            assert(_successful_PRNG_states.size() == 0);
            
        }

    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Path_Policy_Base() {

        }

        /** 
         * record initial state of PRNG if path is successful
         */
        void record_success(const Configuration_type &configuration, const base_generator_type &PRNG_state) {

            if (_successful_PRNG_states.size() < 10) 
                if (success(configuration))
                    _successful_PRNG_states.push_back(PRNG_state);
        }

        /** 
         * return true if path should be terminated
         */
        const bool terminate(const Configuration_type &configuration) const {

            const bool cond1 = configuration.end();
            const bool cond2 = success(configuration);
            const bool cond3 = extinguished(configuration);  // 8th March 2013
            
            return (cond1 || cond2 || cond3);
        }

        /**
         * get collection of initial PRNG states that yield successful trajectories 
         */
        const std::vector<base_generator_type> successful_PRNG_states() const { 
            
            return _successful_PRNG_states;
        }
        
    };

    /**
     * simulate until all sub-populations have extinguished
     */
    template<class Configuration_type>
    class Wait_For_All_Extinctions : public Path_Policy_Base<Configuration_type> {
    private:

        /** 
         * return true if path is successful
         */
        virtual const bool success(const Configuration_type &configuration) const {

            return extinguished(configuration);
            
        }



    };

    /**
     * simulate until last species has arisen in all sub-populations \n
     */
    template<class Configuration_type>
    class Wait_For_Last_Species_In_All_SubPops : public Path_Policy_Base<Configuration_type> {
    private:

        /** 
         * return true if path is successful
         */
        virtual const bool success(const Configuration_type &configuration) const {

            return all_last_species_occurred__within_timeSpan(configuration); // new
            
        }


    };


    /** 
     * terminate path when last species has arisen or sub-pop has extinguished for all sub-pops
     */
    template<class Configuration_type>
    class WaitForEventualFate_LastSpecies_allSubPops : public Path_Policy_Base<Configuration_type> {
    private:

        /** 
         * return true if path is successful
         */
        virtual const bool success(const Configuration_type &configuration) const {

            return eventualFateOccurred_lastSpecies_allSubPops(configuration);
            
        }


    };

    /** 
     * simulate for a fixed duration 
     */
    template<class Configuration_type>
    class Fixed_Duration : public Path_Policy_Base<Configuration_type> {
    private:

        /** 
         * return true if path is successful
         */
        virtual const bool success(const Configuration_type &configuration) const {

            return false;
            
        }


    };

}

#endif	/* PATH_POLICY_H */

