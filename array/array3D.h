#ifndef ARRAY3D_H
#define ARRAY3D_H

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
     * template 3D array\n
     * all members functions are public because this class may be instantiated
     */
    template <class data_type>
    class Array3D {

    protected:

        int dim0, dim1, dim2;
        int dim01;
        std::vector<data_type> data;

        /**
         * check indices \n
         * used only from at() method
         * rely upon std::vector to throw error
         */
        void indices_check(const int &ii, const int &jj, const int &kk) const {

            assert(ii >= 0);
            assert(ii < dim0);
            assert(jj >= 0);
            assert(jj < dim1);
            assert(kk >= 0);
            assert(kk < dim2);
        }


        /**
         * check indices along dim0 \n
         * rely upon std::vector to throw error
         */
        void indices_check0(const int &ii) const {

            assert(ii >= 0);
            assert(ii < dim0);
        }


        /**
         * check indices along dim1 \n
         * rely upon std::vector to throw error
         */
        void indices_check1(const int &jj) const {

            assert(jj >= 0);
            assert(jj < dim1);
        }

        /**
         * check indices along dim3 \n
         * rely upon std::vector to throw error
         */
        void indices_check2(const int &kk) const {

            assert(kk >= 0);
            assert(kk < dim2);
        }

    public:

        /**
         * custom constructor
         */
        explicit Array3D(const int &dim0_, const int &dim1_, const int &dim2_)
        : dim0(dim0_), dim1(dim1_), dim2(dim2_) {

            dim01 = dim0 * dim1;
            data = std::vector<data_type > (dim0 * dim1 * dim2);

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit Array3D() : dim0(0), dim1(0), dim2(0) {

            dim01 = dim0 * dim1;
            data = std::vector<data_type > (dim0 * dim1 * dim2);

        }

        /**
         *  @brief  Provides access to the data contained in a 3D array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read/write reference to data.
         *
         */
        data_type& at(const int &ii, const int &jj, const int &kk) {

#ifdef CHECKBOUNDS
            indices_check(ii, jj, kk);
#endif
            return data.at(ii + dim0 * jj + dim01 * kk);
        }

        /**
         *  @brief  Provides access to the data contained in a <b> const </b> 3D array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read-only reference to data.
         *
         */
        const data_type& at(const int &ii, const int &jj, const int &kk) const {

#ifdef CHECKBOUNDS
            indices_check(ii, jj, kk);
#endif
            return data.at(ii + dim0 * jj + dim01 * kk);
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
         * get number of elements along dimension 2
         */
        const int get_dim2() const {

            return dim2;
        }

        /**
         * get total number of elements
         */
        const int get_dim() const {

            return (dim0 * dim1 * dim2);
        }

    };



} 

#endif
