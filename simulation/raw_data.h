#ifndef RAW_DATA_H
#define	RAW_DATA_H

#include <boost/shared_ptr.hpp> // boost::shared_ptr
#include <boost/type_traits/is_same.hpp>

#include <mutation_accumulation/configuration/configuration/branching_discrete_diamond.h> 

/*************************************************************************/

namespace monte_carlo {

    /** 
     * interface for printing out raw data 
     */
    template <class Configuration_Policy>
    class Raw_Data_Policy_Base {
    protected:

        /**
         * custom constructor
         * this is an empty constructor\n
         * compiler should generate code to construct private data members (Item 30)
         */
        explicit Raw_Data_Policy_Base() {

            /* check template parameter types to supplement "duck typing" */
            typedef typename Configuration_Policy::time_t time_type; // Item 42
            typedef typename Configuration_Policy::population_t population_type;
            BOOST_STATIC_ASSERT((boost::is_base_of<Configuration_Interface<time_type, population_type>, Configuration_Policy>::value));
        }


    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Raw_Data_Policy_Base() {

        }

        /* 
         * print raw data
         */
        virtual void print(const Configuration_Policy &configuration) = 0;

    };

    /** 
     * print out nothing
     */
    template <class Configuration_Policy>
    class Raw_Data_Null : public Raw_Data_Policy_Base<Configuration_Policy> {
    public:

        /* 
         * print raw data
         */
        virtual void print(const Configuration_Policy &configuration) {

        }

    };

    /** 
     * print out first time at which a species-2 stem cell arises in sub-population 0
     */
    template <class Configuration_Policy>
    class Raw_Data_T2 : public Raw_Data_Policy_Base<Configuration_Policy> {
    public:

        const boost::shared_ptr<std::ofstream> _ofstream_ptr;

    public:

        /**
         * custom constructor
         */
        explicit Raw_Data_T2() : _ofstream_ptr(open_file_for_output("T2.dat")) {

        }

        /* 
         * print raw data
         */
        virtual void print(const Configuration_Policy &configuration) {

            *_ofstream_ptr << std::setw(20) << configuration.mutation_times(Pop(0), Spe(2)) << std::endl;

        }

    };

    /** 
     * print out first time at which a double-mutant stem cell arises in sub-population 0 \n
     * also print out whether each single-mutant species generated the first double-mutant \n
     * assumes that Configuration_Policy = Branching_Discrete_Diamond<population_type>, where population_type is variable \n
     */

    template <class Configuration_Policy>
    class Raw_Data_T2_Route : public Raw_Data_Policy_Base<Configuration_Policy > {
    public:

        const boost::shared_ptr<std::ofstream> _ofstream_ptr;

    public:

        /**
         * custom constructor
         */
        explicit Raw_Data_T2_Route() : _ofstream_ptr(open_file_for_output("T2_singleMutant.dat")) {

            /* check template parameter type to supplement duck typing */
            typedef typename Configuration_Policy::population_t population_type;
            BOOST_STATIC_ASSERT((boost::is_same<Branching_Discrete_Diamond<population_type>, Configuration_Policy>::value));

        }

        /* 
         * print raw data
         */
        virtual void print(const Configuration_Policy &configuration) {

            const Spe last_species(configuration.number_species() - 1);                   
            *_ofstream_ptr << std::setw(20) << configuration.mutation_times(Pop(0), last_species);
            *_ofstream_ptr << std::setw(20) << configuration.a_yielded_ab(Pop(0));
            *_ofstream_ptr << std::setw(20) << configuration.b_yielded_ab(Pop(0)) << std::endl;
        }

    };
}



#endif	/* RAW_DATA_H */

