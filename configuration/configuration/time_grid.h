#ifndef TIME_GRID_H
#define	TIME_GRID_H

#include <iomanip> // std::setw 

#include <mutation_accumulation/parameters/parameters_fwd.h> // monte_carlo::Node
#include <mutation_accumulation/utility/grid.h> // grid::make_uniform_grid
#include <mutation_accumulation/simulation/files.h> // monte_carlo::open_file_for_output

#include <boost/current_function.hpp> 

/*************************************************************************/

namespace monte_carlo {

    /**
     * create uniform grid of time points of arbitrary type (int, double) starting at zero\n
     */
    template <class time_type>
    class Uniform_Time_Grid {
    private:

        std::vector<time_type> _node_times;

    private:

        /**
         * print time grid
         */
        void print_time_grid() const {

            /* open file for (over-)writing */
            boost::shared_ptr<std::ofstream> ofstream_ptr = open_file_for_output("time_grid.dat");

            for (int ii = 0; ii < _node_times.size(); ii++)
                *ofstream_ptr << _node_times.at(ii) << std::endl;

        }



    public:

        /**
         * custom constructor\n
         * create time grid using suggested arbitrary number of time points
         */
        explicit Uniform_Time_Grid(const int &suggested_number_time_points, const time_type &time_span) {

            _node_times = grid::make_uniform_grid(suggested_number_time_points, static_cast<time_type>(0), time_span);
            
            print_time_grid();
        }

        /**
         * custom constructor \n
         * create time grid containing just initial and final time points
         */
        explicit Uniform_Time_Grid(const time_type &time_span) {

            _node_times = std::vector<time_type > (2);
            _node_times.at(0) = static_cast<time_type> (0);
            _node_times.at(1) = time_span;
            
            print_time_grid();

        }

        /**
         * default constructor \n
         */
        explicit Uniform_Time_Grid() {

            _node_times = std::vector<time_type > (1);
            _node_times.at(0) = static_cast<time_type> (0);
            
            print_time_grid();

        }

        /**
         * get time at specific node
         */
        const time_type at(const Node &node) const {

            try {
                return _node_times.at(node.value());                
            } catch (const std::out_of_range &ee) { 
                assert(false);
            }
        }

        /**
         * get number of time points
         */
        const int size() const {

            return _node_times.size();
        }

    };


}

#endif	

