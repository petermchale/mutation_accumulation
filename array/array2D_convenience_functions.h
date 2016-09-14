#ifndef ARRAY2D_CONVENIENCE_FUNCTIONS_H
#define ARRAY2D_CONVENIENCE_FUNCTIONS_H

#include <vector> // std::vector

#include "array2D.h"

/*************************************************************************/

namespace array {

    /**
     * sums along the dimension of array_2D specified by int dim\n
     * same interface as Matlab "sum(A, dim)"
     */
    template <class data_type>
    const std::vector<data_type> sum(
    const Array2D<data_type> &array_2D,
    const int &dim) {

        assert(dim >= 0);
        assert(dim < 2);
        
        if (dim == 0) // sum along rows
        {
            std::vector<data_type> sum_along_rows(array_2D.get_dim1(), data_type(0));

            for (int jj = 0; jj < array_2D.get_dim1(); jj++)
                for (int ii = 0; ii < array_2D.get_dim0(); ii++)
                    sum_along_rows.at(jj) += array_2D.at(ii, jj);

            return sum_along_rows;
            
        } else if (dim == 1) // sum along columns
        {
            std::vector<data_type> sum_along_columns(array_2D.get_dim0(), data_type(0));

            for (int ii = 0; ii < array_2D.get_dim0(); ii++)
                for (int jj = 0; jj < array_2D.get_dim1(); jj++)
                    sum_along_columns.at(ii) += array_2D.at(ii, jj);

            return sum_along_columns;
            
        } 
            

    }


} // end namespace array

#endif // ARRAY2D_CONVENIENCE_FUNCTIONS_H
