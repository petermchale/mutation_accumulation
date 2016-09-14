#ifndef POPULATION2D_H
#define	POPULATION2D_H

#include <mutation_accumulation/parameters/parameters_fwd.h> // monte_carlo::Pop, etc
#include <mutation_accumulation/array/array2D.h> // array::Array2D
#include <mutation_accumulation/array/array2D_convenience_functions.h> // array::sum

/*************************************************************************/

namespace monte_carlo {

    /**
     * 2D array to hold population sizes of arbitrary number of species in arbitrary number of sub-populations \n
     * should make Array2D a private data member of Population2D since Population2D doesn't use all public methods of Array2D; Fowler p352
     * 
     */
    template<class population_type>
    class Population2D : public array::Array2D<population_type> {
    public:

        typedef population_type population_t;

    private:

        static const array::Array2D<population_type> vector_to_Array2D(const std::vector<population_type> &input_vector) {

            const int number_rows = 1;
            const int number_cols = input_vector.size();
            array::Array2D<population_type> array2D(number_rows, number_cols);
            for (int col = 0; col < number_cols; col++)
                array2D.at(0, col) = input_vector.at(col);

            return array2D;
        }


    public:

        /**
         * custom constructor
         */
        explicit Population2D(const array::Array2D<population_type> &array2D)
        : array::Array2D<population_type>(array2D) {

        }

        /**
         * custom constructor
         */
        explicit Population2D(const int &number_sub_pops_, const int &number_species_)
        : array::Array2D<population_type>(number_sub_pops_, number_species_) {

        }

        /**
         * custom constructor that takes vector of species in a single subpopulation as argument
         */
        explicit Population2D(const std::vector<population_type> &input_vector)
        : array::Array2D<population_type>(vector_to_Array2D(input_vector)) {

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit Population2D() : array::Array2D<population_type>() {

        }

        /**
         * number of sub-populations
         */
        const int number_sub_pops() const {

            return (this->dim0);
        }

        /**
         * number of species
         */
        const int number_species() const {

            return (this->dim1);
        }

    };

    /**
     * calculate size of sub-population pop in a Population2D object
     */
    template <class population_type>
    const population_type sub_pop_size(
    const Population2D<population_type> &populations,
    const Pop &pop) {

        const int column_dim = 1;

        // sum along columns using row defined by pop
        return array::sum(populations, column_dim).at(pop.value());

    }



}

#endif	/* POPULATION2D_H */

