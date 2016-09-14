#ifndef ARRAY4D_CONVENIENCE_FUNCTIONS_H
#define ARRAY4D_CONVENIENCE_FUNCTIONS_H

#include "array4D.h"
#include <fstream> // std::ofstream

/*************************************************************************/

namespace array {

    /**
     * store 4D array to disk for inspection
     */
    template <class data_type>
    void store_4D_array(
    const Array4D<data_type> &array_4D,
    const std::string &fileName) {

        std::ofstream fout(fileName.c_str()); // file closed when fout goes out of scope (RAII)

        if (!fout) {
            std::cerr << "cannot open " << fileName << std::endl;
            assert(fout);
        }

        /* executing the loops in this order ensures that
         * we access contiguous parts of memory */
        for (int ll = 0; ll < array_4D.get_dim3(); ll++)
            for (int kk = 0; kk < array_4D.get_dim2(); kk++)
                for (int jj = 0; jj < array_4D.get_dim1(); jj++)
                    for (int ii = 0; ii < array_4D.get_dim0(); ii++)
                        fout << std::setw(10) << std::setprecision(3) << array_4D.at(ii, jj, kk, ll);

        fout << std::endl;
    }

} // end namespace array

#endif // ARRAY4D_CONVENIENCE_FUNCTIONS_H
