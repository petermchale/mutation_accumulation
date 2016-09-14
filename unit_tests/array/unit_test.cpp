#include <iostream> // std::cout
#include <iomanip> // std::setw
#include <cassert> // assert

#include <mutation_accumulation/array/array2D.h> // array:Array2D
#include <mutation_accumulation/array/array2D_convenience_functions.h> // array::sum

#include "unit_test.h"

/*************************************************************************/

void unit_test::unit_test_array_sum_2D() {

    std::cout << "testing array::sum on 2D array..." << std::endl;
    
    array::Array2D<int> array_2D(2, 3);

    /* load 2D array */
    array_2D.at(0, 0) = 1;
    array_2D.at(0, 1) = 2;
    array_2D.at(0, 2) = 3;

    array_2D.at(1, 0) = 4;
    array_2D.at(1, 1) = 5;
    array_2D.at(1, 2) = 6;

//    /* print 2D array in human-readable format */
//    const int column_width = 5;
//    std::cout << "array = " << std::endl;
//    std::cout << std::setw(column_width) << array_2D.at(0, 0);
//    std::cout << std::setw(column_width) << array_2D.at(0, 1);
//    std::cout << std::setw(column_width) << array_2D.at(0, 2);
//    std::cout << std::endl;
//    std::cout << std::setw(column_width) << array_2D.at(1, 0);
//    std::cout << std::setw(column_width) << array_2D.at(1, 1);
//    std::cout << std::setw(column_width) << array_2D.at(1, 2);
//    std::cout << std::endl;

    {
        /* check summing along rows */

        std::cout << "checking summing along rows..." << std::endl;
        const int dim = 0;
        assert(array::sum(array_2D, dim).size() == array_2D.get_dim1());
        std::cout << "passed: function returns vector whose size is the number of columns in 2D array" << std::endl;

        assert(array::sum(array_2D, dim).at(0) == (array_2D.at(0, 0) + array_2D.at(1, 0)));
        assert(array::sum(array_2D, dim).at(1) == (array_2D.at(0, 1) + array_2D.at(1, 1)));
        assert(array::sum(array_2D, dim).at(2) == (array_2D.at(0, 2) + array_2D.at(1, 2)));
        std::cout << "passed: function correctly sums along rows for each column" << std::endl;
    }

    {
        /* check summing along columns */

        std::cout << "checking summing along columns..." << std::endl;
        const int dim = 1;
        assert(array::sum(array_2D, dim).size() == array_2D.get_dim0());
        std::cout << "passed: function returns vector whose size is the number of rows in 2D array" << std::endl;

        assert(array::sum(array_2D, dim).at(0) == (array_2D.at(0, 0) + array_2D.at(0, 1) + array_2D.at(0,2)));
        assert(array::sum(array_2D, dim).at(1) == (array_2D.at(1, 0) + array_2D.at(1, 1) + array_2D.at(1,2)));
        std::cout << "passed: function correctly sums along columns for each row" << std::endl;
    }


    std::cout << "...all tests passed" << std::endl;
    
}


