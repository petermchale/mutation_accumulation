#ifndef PATH_H
#define	PATH_H

#include "population2D.h" // monte_carlo::Population2D
#include "time_grid.h" // Uniform_Time_Grid

/*************************************************************************/

namespace monte_carlo {

    /**
     * store configurations at a grid of time points 
     */
    template <class time_type, class population_type>
    class Path {
    private:

        typedef Population2D<population_type> population2D_type;
        std::vector<population2D_type> _node_values;
        Uniform_Time_Grid<time_type> _node_times; // must contain at least the initial time point and the final time point

    public:

        /**
         * custom constructor
         */
        explicit Path(
                const population2D_type &initial_node_value,
                const Uniform_Time_Grid<time_type> &node_times)
        :
        _node_values(std::vector<population2D_type>(1, initial_node_value)),
        _node_times(node_times) {

        }

        /**
         * is path complete?
         */
        const bool complete() const {

            const int number_filled_nodes = _node_values.size();
            const int total_number_nodes = _node_times.size();

            if (number_filled_nodes < total_number_nodes) {
                return false;
            } else if (number_filled_nodes == total_number_nodes) {
                return true;
            } else {
                std::cerr << "number_filled_nodes > total_number_nodes" << std::endl;
                assert(false);
            }

        }

        /**
         * update path if current time has reached next node \n
         * this implementation valid for discrete time only \n
         */
        void update(const int &tt, const population2D_type &value) {

            const Node next_node(_node_values.size());

            if (tt >= _node_times.at(next_node))
                _node_values.push_back(value);

        }

        /**
         * update path \n
         * this implementation valid for continuous time only \n
         */
        void update(const double &tt, const population2D_type &value, const population2D_type &value_old) {

            while (!complete()) {
                const Node next_node(_node_values.size());
                if (tt > _node_times.at(next_node))
                    _node_values.push_back(value_old);
                else
                    return;
            }
        }

        /**
         * get population of species spe in sub-population pop at time indicated by node
         */
        const population_type at(const Pop &pop, const Spe &spe, const Node &node) const {

            population2D_type node_value;
            try {
                node_value = _node_values.at(node.value());
            } catch (const std::out_of_range &ee) {
                assert(false);
            }

            return node_value.at(pop.value(), spe.value());
        }

        /**
         * get time indicated by node
         */
        const time_type get_time(const Node &node) const {

            return _node_times.at(node);
        }

        /**
         * number of filled nodes 
         */
        const int number_filled_nodes() const {

            return _node_values.size();
        }

        /**
         * time value of last node
         */
        const time_type last_node_time() const {

            const Node last_node(_node_times.size() - 1);

            return _node_times.at(last_node);

        }



    };



}

#endif	/* PATH_H */

