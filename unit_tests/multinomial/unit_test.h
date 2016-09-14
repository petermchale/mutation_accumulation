#ifndef UNIT_TEST_H
#define	UNIT_TEST_H

/*************************************************************************/

namespace unit_test {

    /**
     * 2D multinomial should coincide with binomial CDF from matlab\n
     */
    void unit_test_2Dmultinomial();

    /**
     * 2D multinomial should produce  a vector whose elements sume to N\n
     */
    void unit_test__sum_random_vector_elements();


    /**
     * 3D multinomial should coincide with trinomial PDF from matlab\n
     */
    void unit_test__3Dmultinomial();


    /**
     * 4D multinomial should coincide with PDF from matlab\n
     */
    void unit_test__4Dmultinomial();



    /**
     * 3D multinomial should always generate a single random vector when all but one of the categorical probabilities are zero\n
     */
    void unit_test__3Dmultinomial_edgeCase();


    /**
     * 3D multinomial should always generate (0,0,0) when N = 0\n
     */
    void unit_test__3Dmultinomial_zeroN();


}

#endif	/* UNIT_TEST_H */

