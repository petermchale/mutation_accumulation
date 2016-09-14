#ifndef ARRAY2D_H
#define ARRAY2D_H

#define CHECKBOUNDS

#include <assert.h> // assert
#include <deque> // std::deque

/**
 * Code is based upon Philippe Monthoux's and Clare Yu's \n
 * Also modelled on stl_vector.h
 * 
 * One could derive eg 2D, 3D, 4D array classes from a ND array base class \n
 * that stores the size along each dimension in a std::container of size N\n
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
     * template 2D array\n
     * all members functions are public because this class may be instantiated \n
     * uses deque to handle bool template parameter \n
     * see : http://stackoverflow.com/questions/3458856/overloaded-operator-on-template-class-in-c-with-const-nonconst-versions
     */
    template <class data_type>
    class Array2D {
    protected:

        int dim0, dim1;
        std::deque<data_type> data;

        /**
         * check indices \n
         * used only from at() method
         * rely upon std::deque to throw error
         */
        void indices_check(const int &ii, const int &jj) const {

            assert(ii >= 0);
            assert(ii < dim0);
            assert(jj >= 0);
            assert(jj < dim1);
        }

        /**
         * check indices along dim0 \n
         * rely upon std::deque to throw error
         */
        void indices_check0(const int &ii) const {

            assert(ii >= 0);
            assert(ii < dim0);
        }

        /**
         * check indices along dim1 \n
         * rely upon std::deque to throw error
         */
        void indices_check1(const int &jj) const {

            assert(jj >= 0);
            assert(jj < dim1);
        }

    public:

        /**
         * custom constructor
         */
        explicit Array2D(const int &dim0_, const int &dim1_)
        : dim0(dim0_), dim1(dim1_) {

            data = std::deque<data_type > (dim0 * dim1);

        }

        /**
         * custom constructor
         */
        explicit Array2D(const int &dim0_, const int &dim1_, const data_type &initial_value)
        : dim0(dim0_), dim1(dim1_) {

            data = std::deque<data_type > (dim0 * dim1, initial_value);

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit Array2D() : dim0(0), dim1(0) {

            data = std::deque<data_type > (dim0 * dim1);

        }

        /**
         *  @brief  Provides access to the data contained in a 2D array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read/write reference to data.
         *
         */
        data_type& at(const int &ii, const int &jj) {

#ifdef CHECKBOUNDS
            indices_check(ii, jj);
#endif
            return data.at(ii + dim0 * jj);
        }

        /**
         *  @brief  Provides access to the data contained in a <b> const </b> 2D array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read-only reference to data.
         *
         */
        const data_type& at(const int &ii, const int &jj) const {

#ifdef CHECKBOUNDS
            indices_check(ii, jj);
#endif
            return data.at(ii + dim0 * jj);
        }


        /**
         * get number of elements along dimension 0
         */
        const int get_dim0() const {

            return dim0;
        }

        /**
         * get number of elements along dimension 1
         */
        const int get_dim1() const {

            return dim1;
        }

        /**
         * get total number of elements
         */
        const int get_dim() const {

            return (dim0 * dim1);
        }

        /**
         * calculate sum of all elements in 2D array\n
         *
         * kept this as member function for the sake of performance
         */
        const data_type sum() const {

            data_type sum = data_type(0);

            for (int nn = 0; nn < data.size(); nn++)
                sum += data.at(nn);

            return sum;
        }


    };



} 

#endif
