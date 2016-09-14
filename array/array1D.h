#ifndef ARRAY1D_H
#define ARRAY1D_H

#define CHECKBOUNDS

#include <assert.h> // assert
#include <vector> // std::vector

/**
 * Code is based upon Philippe Monthoux's and Clare Yu's \n
 * Also modelled on stl_vector.h
 * 
 * One could derive eg 2D, 3D, 4D array classes from a ND array base class \n
 * that stores the size along each dimension in a std::vector of size N\n
 * sum method could then be implemented in the base class
 *
 * Other multi-dimensional array libraries to consider include \n
 * Boost: but overly complex syntax
 * Blitz++: supposed to be fast\n
 * lite array: nice syntax but seems to be used by no-one and un-maintained (at least since 2009)
 */

/*************************************************************************/

/**
 * this namespace includes classes that store multi-dimensional arrays\n
 * should really make protected data members private (Item 22)
 */
namespace array {

    /**
     * template 1D array\n
     * all members functions are public because this class may be instantiated
     */
    template <class data_type>
    class Array1D {

    protected:

        int dim;
        std::vector<data_type> data;

        /**
         * check indices \n
         * used only from at() method
         * rely upon std::vector to throw error
         */
        void indices_check(const int &ii) const {

            assert(ii >= 0);
            assert(ii < dim);
        }

    public:

        /**
         * custom constructor
         */
        explicit Array1D(const int &dim_)
        : dim(dim_) {

//            data = std::vector<data_type > (dim, data_type(0));
            data = std::vector<data_type > (dim);

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit Array1D() : dim(0) {

//            data = std::vector<data_type > (dim, data_type(0));
            data = std::vector<data_type > (dim);

        }

        /**
         *  @brief  Provides access to the data contained in an array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read/write reference to data.
         *
         */
        data_type& at(const int &ii) {

#ifdef CHECKBOUNDS
            indices_check(ii);
#endif
            return data.at(ii);
        }

        /**
         *  @brief  Provides access to the data contained in a <b> const </b> array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read-only reference to data.
         *
         */
        const data_type& at(const int &ii) const {

#ifdef CHECKBOUNDS
            indices_check(ii);
#endif
            return data.at(ii);
        }

        /**
         * get total number of elements
         */
        const int get_dim() const {

            return dim;
        }

    };



} // end namespace array

#endif
