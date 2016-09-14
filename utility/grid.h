#ifndef GRID_H
#define	GRID_H

#include <cassert> // assert
#include <vector> // std::vector
#include <algorithm>    // std::sort; std::unique
#include <cmath> // pow
#include <iostream> // cerr 

#include <mutation_accumulation/utility/data_traits.h> // data_types::discrete_type

/*************************************************************************/

namespace grid {

    namespace detail {

        /** 
         * create uniform grid vector of arbitrary type 
         */
        template <class value_type>
        const std::vector<value_type> build_grid_vector(
        const int &grid_size,
        const value_type &delta_grid,
        const value_type &grid_lower) {

            std::vector<value_type> grid(grid_size, static_cast<value_type> (-1));
            for (int ii = 0; ii < grid.size(); ii++)
                grid.at(ii) = grid_lower + static_cast<value_type> (ii) * delta_grid;

            return grid;
        }

        /**
         * create uniform grid of discrete type (int or long long int)\n
         */
        template <class value_type>
        const std::vector<value_type> make_uniform_grid(
        const int &suggested_number_grid_points,
        const value_type &grid_lower,
        const value_type &grid_upper,
        data_types::discrete_type) {

            /* size of return vector */
            int vector_size = suggested_number_grid_points;
            int vector_size_bar = suggested_number_grid_points - 1;

            /* grid span */
            const value_type span_of_grid = grid_upper - grid_lower;

            /* size of increment in vector elements */
            value_type delta_grid = -1;
            if (vector_size_bar <= span_of_grid) {

                if ((span_of_grid % vector_size_bar) == 0) {

                    delta_grid = span_of_grid / vector_size_bar;

                } else {

                    while ((span_of_grid % vector_size_bar) != 0) {
                        vector_size++;
                        vector_size_bar++;
                    }

                    delta_grid = span_of_grid / vector_size_bar;
                }

            } else {

                vector_size = span_of_grid + 1;
                vector_size_bar = span_of_grid;
                delta_grid = 1;
            }

            /* use vector size and vector increment to build return vector */
            return build_grid_vector(vector_size, delta_grid, grid_lower);

        }

        /**
         * create uniform grid of continuous type (double) \n
         */
        template <class value_type>
        const std::vector<value_type> make_uniform_grid(
        const int &number_grid_points,
        const value_type &grid_lower,
        const value_type &grid_upper,
        data_types::continuous_type) {

            /* size of return vector */
            int vector_size = number_grid_points;
            int vector_size_bar = number_grid_points - 1;

            /* distance between grid points */
            const value_type delta_grid = (grid_upper - grid_lower) / static_cast<value_type> (vector_size_bar);

            /* use vector size and vector increment to build return vector */
            return build_grid_vector(vector_size, delta_grid, grid_lower);
        }

    }

    namespace detail {

        /** 
         * remove duplicate entries from vector\n
         * works for entries of discrete type
         */
        template <class value_type>
        const std::vector<value_type> remove_duplicate_elements(const std::vector<value_type> &data_in) {

            std::vector<value_type> data_out = data_in;
            std::sort(data_out.begin(), data_out.end());
            data_out.erase(std::unique(data_out.begin(), data_out.end()), data_out.end());

            return data_out;
        }

        /**
         * create logarithmic grid of arbitrary type (int, long long int, double)\n
         * may create duplicate elements if value_type is discrete \n
         */
        template <class value_type>
        const std::vector<value_type> do_make_logarithmic_grid(
        const int &number_grid_points,
        const value_type &grid_lower,
        const value_type &grid_upper) {

            if (grid_lower > static_cast<value_type> (0)) {

                /* create cts uniform grid in log space */
                const std::vector<double> grid_exponents = detail::make_uniform_grid(number_grid_points, log10(grid_lower), log10(grid_upper), data_types::continuous_type());

                /* convert back to original space and original type */
                std::vector<value_type> grid(grid_exponents.size());
                for (int ii = 0; ii < grid.size(); ii++)
                    grid.at(ii) = static_cast<value_type> (pow(10, grid_exponents.at(ii)));

                return grid;

            } else {

                std::cerr << "lower bound of sample space must be greater than zero if log is to be taken" << std::endl;
                assert(false);
            }


        }

        /**
         * create logarithmic grid of discrete type (int, long long int)\n
         */
        template <class value_type>
        const std::vector<value_type> make_logarithmic_grid(
        const int &number_grid_points,
        const value_type &grid_lower,
        const value_type &grid_upper,
        data_types::discrete_type) {

            return remove_duplicate_elements(do_make_logarithmic_grid(number_grid_points, grid_lower, grid_upper));
            
        }

        /**
         * create logarithmic grid of continuous type (double)\n
         */
        template <class value_type>
        const std::vector<value_type> make_logarithmic_grid(
        const int &number_grid_points,
        const value_type &grid_lower,
        const value_type &grid_upper,
        data_types::continuous_type) {

            return do_make_logarithmic_grid(number_grid_points, grid_lower, grid_upper);

        }


    }

    /**
     * create uniform grid of arbitrary type (int, long long int, double)\n
     */
    template <class value_type>
    const std::vector<value_type> make_uniform_grid(
    const int &suggested_number_grid_points,
    const value_type &grid_lower,
    const value_type &grid_upper) {

        typename data_types::data_traits<value_type>::category value_category;
        return detail::make_uniform_grid(suggested_number_grid_points, grid_lower, grid_upper, value_category);
    }

    /**
     * create logarithmic grid of arbitrary type (int, long long int, double)\n
     */
    template <class value_type>
    const std::vector<value_type> make_logarithmic_grid(
    const int &number_grid_points,
    const value_type &grid_lower,
    const value_type &grid_upper) {

        typename data_types::data_traits<value_type>::category value_category;
        return detail::make_logarithmic_grid(number_grid_points, grid_lower, grid_upper, value_category);

    }


}

#endif	/* GRID_H */


