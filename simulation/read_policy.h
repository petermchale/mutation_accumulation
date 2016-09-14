#ifndef READ_POLICY_H
#define	READ_POLICY_H

//#define DEBUG_READ_POLICY_H

#include <boost/shared_ptr.hpp> // boost::shared_ptr
#include <mutation_accumulation/parameters/parameters_fwd.h> // Symmetry, etc

#include <boost/type_traits/is_base_of.hpp>
#include <mutation_accumulation/configuration/configuration/configuration_interface.h>
#include <mutation_accumulation/configuration/configuration/population2D.h>
#include <mutation_accumulation/configuration/configuration/mutation_rates.h> // monte_carlo::MutationRates
#include <mutation_accumulation/configuration/configuration/time_grid.h> // monte_carlo::Uniform_Time_Grid

#include <mutation_accumulation/utility/strings.h> // strings::parse_scalar, etc

/*************************************************************************/

/** 
 * read in parameter values from input file 
 */
namespace monte_carlo {

    /** 
     * interface for reading in parameter values
     */
    template <class Configuration_Policy>
    class Read_Policy_Base {
    public:
        typedef typename Configuration_Policy::time_t time_t; // Item 42
        typedef typename Configuration_Policy::population_t population_t;
        typedef Population2D<population_t> Population2D_t;
        typedef Uniform_Time_Grid<time_t> Uniform_Time_Grid_t;
        typedef double error_t;
        typedef int divisor_t;

    protected:

        /**
         * custom constructor
         * this is an empty constructor\n
         * compiler should generate code to construct private data members (Item 30)
         */
        explicit Read_Policy_Base() {

            /* check template parameter types to supplement "duck typing" */
            BOOST_STATIC_ASSERT((boost::is_base_of<Configuration_Interface<time_t, population_t>, Configuration_Policy>::value));
        }


    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Read_Policy_Base() {

        }

        /**
         * get population size
         */
        virtual const Population2D_t get_population() const = 0;

        /**
         * get renewal bias
         */
        virtual const SymmetricRenewal get_SymmetricRenewal() const = 0;

        /**
         * get mutation rate
         */
        virtual const MutationRates getUu() const = 0;

        /**
         * get symmetry
         */
        virtual const Symmetry getSymmetry() const = 0;

        /**
         * get time grid
         */
        virtual const Uniform_Time_Grid_t getTime_grid() const = 0;

        /**
         * get time span for histogram
         */
        virtual const time_t getTime_span_histogram() const = 0;

        /**
         * get error in probability 
         */
        virtual const error_t get_error_prob() const = 0;

        /**
         * get observer divisor 
         */
        virtual const divisor_t get_observer_divisor() const = 0;

    };

    /** 
     * read in parameters, assuming renewal bias = 0.5 (homeostasis)
     */
    template <class Configuration_Policy>
    class Read_Homeostasis_Policy : public Read_Policy_Base<Configuration_Policy > {
    private:
        typedef typename Configuration_Policy::time_t time_type; // Item 42
        typedef typename Configuration_Policy::population_t population_type;
        typedef Population2D<population_type> Population2D_type;
        typedef Uniform_Time_Grid<time_type> Uniform_Time_Grid_type;
        typedef typename Read_Policy_Base<Configuration_Policy>::error_t error_type; // item 43 "Know how to access names in templatized base classes" in "Effective C++: 55 Specific Ways to Improve Your Programs and Designs"
        typedef typename Read_Policy_Base<Configuration_Policy>::divisor_t divisor_type; // item 43 "Know how to access names in templatized base classes" in "Effective C++: 55 Specific Ways to Improve Your Programs and Designs"

        Population2D_type _population2D;
        MutationRates _uu;
        Symmetry _symmetry;
        Uniform_Time_Grid_type _time_grid;
        time_type _time_span_histogram;
        error_type _error_probability;
        divisor_type _observer_divisor;


    public:

        /**
         * custom constructor
         */
        explicit Read_Homeostasis_Policy(const std::string &filename) {

            /* open input file */
            const boost::shared_ptr<std::ifstream> ifstream_ptr = open_file_for_input(filename);

            /* read in initial population sizes */
            _population2D = Population2D_type(strings::parse_matrix<population_type> (*ifstream_ptr));

#ifdef DEBUG_READ_POLICY_H
            /* check _population2D */
            for (int pop = 0; pop < _population2D.number_sub_pops(); pop++) {
                for (int spe = 0; spe < _population2D.number_species(); spe++)
                    std::cout << "<" << _population2D.at(pop, spe) << ">" << " ";
                std::cout << std::endl;
            }
#endif

            /* read in mutation rates */
            _uu = MutationRates(strings::parse_vector<MutationRates::data_t > (*ifstream_ptr));

#ifdef DEBUG_READ_POLICY_H 
            /* check _uu */
            for (int ii = 0; ii < _uu.get_dim(); ii++)
                std::cout << "<" << _uu.at(ii) << ">" << std::endl;
#endif

            /* read in symmetry */
            _symmetry = Symmetry(strings::parse_scalar<Symmetry > (*ifstream_ptr));

            /* read in time span for path */
            const time_type time_span_path = strings::parse_scalar<time_type > (*ifstream_ptr);

            /* create a uniform grid of time points at which to sample configuration */
            // Uniform_Time_Grid<time_type> time_grid(suggested_number_time_points, time_span_path);
            _time_grid = Uniform_Time_Grid_type(time_span_path);

            /* read in span of sample space for histogram */
            _time_span_histogram = strings::parse_scalar<time_type > (*ifstream_ptr);

            /* read in error in probability */
            _error_probability = strings::parse_scalar<error_type> (*ifstream_ptr);

            /* write histogram log data every observer_divisor times observer is notified */
            _observer_divisor = strings::parse_scalar<divisor_type > (*ifstream_ptr);


        }

        /**
         * get population size
         */
        virtual const Population2D_type get_population() const {

            return _population2D;
        }

        /**
         * get renewal bias
         */
        virtual const SymmetricRenewal get_SymmetricRenewal() const {

            return SymmetricRenewal(0.5);
        }

        /**
         * get mutation rate
         */
        virtual const MutationRates getUu() const {

            return _uu;
        }

        /**
         * get symmetry
         */
        virtual const Symmetry getSymmetry() const {
            return _symmetry;
        }

        /**
         * get time grid
         */
        virtual const Uniform_Time_Grid_type getTime_grid() const {
            return _time_grid;
        }

        /**
         * get time span for histogram
         */
        virtual const time_type getTime_span_histogram() const {
            return _time_span_histogram;
        }

        /**
         * get error in probability 
         */
        virtual const error_type get_error_prob() const {
            return _error_probability;
        }

        /**
         * get observer divisor 
         */
        virtual const divisor_type get_observer_divisor() const {
            return _observer_divisor;
        }

    };

    /** 
     * read in parameters, arbitrary renewal bias (e.g. growth)
     */
    template <class Configuration_Policy>
    class Read_NonHomeostasis_Policy : public Read_Policy_Base<Configuration_Policy > {
    private:
        typedef typename Configuration_Policy::time_t time_type; // Item 42
        typedef typename Configuration_Policy::population_t population_type;
        typedef Population2D<population_type> Population2D_type;
        typedef Uniform_Time_Grid<time_type> Uniform_Time_Grid_type;
        typedef typename Read_Policy_Base<Configuration_Policy>::error_t error_type; // item 43 "Know how to access names in templatized base classes" in "Effective C++: 55 Specific Ways to Improve Your Programs and Designs"
        typedef typename Read_Policy_Base<Configuration_Policy>::divisor_t divisor_type; // item 43 "Know how to access names in templatized base classes" in "Effective C++: 55 Specific Ways to Improve Your Programs and Designs"

        Population2D_type _population2D;
        SymmetricRenewal _symmetric_renewal;
        MutationRates _uu;
        Symmetry _symmetry;
        Uniform_Time_Grid_type _time_grid;
        time_type _time_span_histogram;
        error_type _error_probability;
        divisor_type _observer_divisor;
        


    public:

        /**
         * custom constructor
         */
        explicit Read_NonHomeostasis_Policy(const std::string &filename) {

            /**
             * non-homeostatic case is implemented only as a discrete-time branching process (so far) 
             * check template parameter type to supplement duck typing 
             */
            BOOST_STATIC_ASSERT((boost::is_same<Branching_Discrete<population_type>, Configuration_Policy>::value));

            /* open input file */
            const boost::shared_ptr<std::ifstream> ifstream_ptr = open_file_for_input(filename);

            /* read in initial population sizes */
            _population2D = Population2D_type(strings::parse_matrix<population_type> (*ifstream_ptr));

#ifdef DEBUG_READ_POLICY_H
            /* check _population2D */
            for (int pop = 0; pop < _population2D.number_sub_pops(); pop++) {
                for (int spe = 0; spe < _population2D.number_species(); spe++)
                    std::cout << "<" << _population2D.at(pop, spe) << ">" << " ";
                std::cout << std::endl;
            }
#endif

            /* read in symmetric renewal */
            _symmetric_renewal = SymmetricRenewal(strings::parse_scalar<SymmetricRenewal > (*ifstream_ptr));

            /* read in mutation rates */
            _uu = MutationRates(strings::parse_vector<MutationRates::data_t > (*ifstream_ptr));

#ifdef DEBUG_READ_POLICY_H 
            /* check _uu */
            for (int ii = 0; ii < _uu.get_dim(); ii++)
                std::cout << "<" << _uu.at(ii) << ">" << std::endl;
#endif

            /* read in symmetry */
            _symmetry = Symmetry(strings::parse_scalar<Symmetry > (*ifstream_ptr));

            /* read in time span for path */
            const time_type time_span_path = strings::parse_scalar<time_type > (*ifstream_ptr);

            /* create a uniform grid of time points at which to sample configuration */
            // Uniform_Time_Grid<time_type> time_grid(suggested_number_time_points, time_span_path);
            _time_grid = Uniform_Time_Grid_type(time_span_path);

            /* read in span of sample space for histogram */
            _time_span_histogram = strings::parse_scalar<time_type > (*ifstream_ptr);

            /* read in error in probability */
            _error_probability = strings::parse_scalar<error_type> (*ifstream_ptr);

            /* write histogram log data every observer_divisor times observer is notified */
            _observer_divisor = strings::parse_scalar<divisor_type > (*ifstream_ptr);


        }

        /**
         * get population size
         */
        virtual const Population2D_type get_population() const {

            return _population2D;
        }

        /**
         * get renewal bias
         */
        virtual const SymmetricRenewal get_SymmetricRenewal() const {

            return _symmetric_renewal;
        }

        /**
         * get mutation rate
         */
        virtual const MutationRates getUu() const {

            return _uu;
        }

        /**
         * get symmetry
         */
        virtual const Symmetry getSymmetry() const {
            return _symmetry;
        }

        /**
         * get time grid
         */
        virtual const Uniform_Time_Grid_type getTime_grid() const {
            return _time_grid;
        }

        /**
         * get time span for histogram
         */
        virtual const time_type getTime_span_histogram() const {
            return _time_span_histogram;
        }

        /**
         * get error in probability 
         */
        virtual const error_type get_error_prob() const {
            return _error_probability;
        }

        /**
         * get observer divisor 
         */
        virtual const divisor_type get_observer_divisor() const {
            return _observer_divisor;
        }

    };
}

#endif	/* READ_POLICY_H */

