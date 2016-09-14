#ifndef ARRAY4D_H
#define ARRAY4D_H

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
     * template 4D array\n
     * all members functions are public because this class may be instantiated
     */
    template <class data_type>
    class Array4D {
    protected:

        int dim0, dim1, dim2, dim3;
        int dim01, dim012;
        std::vector<data_type> data;

        /**
         * check indices \n
         * used only from at() method
         * rely upon std::vector to throw error
         */
        void indices_check(const int &ii, const int &jj, const int &kk, const int &ll) const {

            assert(ii >= 0);
            assert(ii < dim0);
            assert(jj >= 0);
            assert(jj < dim1);
            assert(kk >= 0);
            assert(kk < dim2);
            assert(ll >= 0);
            assert(ll < dim3);

        }

    public:

        /**
         * custom constructor
         */
        explicit Array4D(const int &dim0_, const int &dim1_, const int &dim2_, const int &dim3_)
        : dim0(dim0_), dim1(dim1_), dim2(dim2_), dim3(dim3_) {

            dim01 = dim0 * dim1;
            dim012 = dim0 * dim1 * dim2;
//            data = std::vector<data_type > (dim0 * dim1 * dim2 * dim3, data_type(0));
            data = std::vector<data_type > (dim0 * dim1 * dim2 * dim3);

        }

        /**
         * must declare default constructor: \n
         * compiler provides default constructor only if no custom constructors exist.
         */
        explicit Array4D() : dim0(0), dim1(0), dim2(0), dim3(0) {

            dim01 = dim0 * dim1;
            dim012 = dim0 * dim1 * dim2;
//            data = std::vector<data_type > (dim0 * dim1 * dim2 * dim3, data_type(0));
            data = std::vector<data_type > (dim0 * dim1 * dim2 * dim3);

        }

        /**
         *  @brief  Provides access to the data contained in a 4D array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read/write reference to data.
         *
         */
        data_type& at(const int &ii, const int &jj, const int &kk, const int &ll) {

#ifdef CHECKBOUNDS
            indices_check(ii, jj, kk, ll);
#endif

            return data.at(ii + dim0 * jj + dim01 * kk + dim012 * ll);
        }

        /**
         *  @brief  Provides access to the data contained in a <b> const </b> 4D array.
         *  @param  The indices of the element for which data should be accessed.
         *  @return  Read-only reference to data.
         *
         */
        const data_type& at(const int &ii, const int &jj, const int &kk, const int &ll) const {

#ifdef CHECKBOUNDS
            indices_check(ii, jj, kk, ll);
#endif

            return data.at(ii + dim0 * jj + dim01 * kk + dim012 * ll);
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
         * get number of elements along dimension 3
         */
        const int get_dim3() const {

            return dim3;
        }

        /**
         * get total number of elements in array
         */
        const int get_dim() const {

            return (dim0 * dim1 * dim2 * dim3);
        }

        /**
         * calculate sum of all elements in 4D array\n
         *
         * Decided to keep this as a member function \n
         * (as opposed to a non-member function)\n
         * because it reduces the number of arithmetical operations\n
         * and makes code faster
         */
        const data_type sum() const {

            data_type sum = data_type(0);

            for (int nn = 0; nn < data.size(); nn++)
                sum += data.at(nn);

            return sum;
        }

        /**
         * find indices at which cumulative sum first exceeds fraction_sum \n
         *
         * function returns true if fraction_sum is less than
         * sum of all elements in the array
         * and false otherwise\n
         *
         * Decided to keep this as a member function \n
         * (as opposed to a non-member function)\n
         * because it reduces the number of arithmetical operations\n
         * and makes code faster
         */
        const bool cumulative_sum(
                const data_type &fraction_sum,
                int &i_cum, int &j_cum, int &k_cum, int &l_cum) const {

            /**
             * set indices to nonsense values\n
             * indices remain set to nonsense values if function returns false
             */
            i_cum = -1;
            j_cum = -1;
            k_cum = -1;
            l_cum = -1;

            data_type running_fraction_sum = data_type(0);

            int nn = 0;

            /* executing the for loops in this order ensures that
             * we access contiguous parts of memory */
            for (int ll = 0; ll < dim3; ll++)
                for (int kk = 0; kk < dim2; kk++)
                    for (int jj = 0; jj < dim1; jj++)
                        for (int ii = 0; ii < dim0; ii++, nn++) {

                            running_fraction_sum += data.at(nn);

                            if (fraction_sum < running_fraction_sum) {

                                i_cum = ii;
                                j_cum = jj;
                                k_cum = kk;
                                l_cum = ll;

                                return true; // this breaks multiple nested for loops once indices are found
                            }
                        }

            return false; // this statement is reached if fraction_sum > sum of all elements (assuming sum of all elements > 0)

        }

    };

} // end namespace array

#endif
