#include <cassert> // assert

#include <boost/random/mersenne_twister.hpp> // boost::random::mt19937
#include <boost/assign/list_of.hpp> // boost::assign::list_of()

#include <mutation_accumulation/random/multinomial_distribution.h> // random::multinomial_distribution
#include <mutation_accumulation/probability/cdf.h> // probability::CDF_discrete

#include "unit_test.h"

typedef int int_type; // probability::CDF_discrete complains at compile-time if int_type is not int (for example long long int)
typedef mutation_accumulation::random::multinomial_distribution<int_type, double> mn_type;

/*************************************************************************/


void unit_test::unit_test_2Dmultinomial() {

    std::cout << "generate 2D multinomial CDF..." << std::endl;

    /* PRNG engine */
    boost::random::mt19937 urng;

    /* parameters of multinomial_distribution */
    const int_type NN = 10;
    const double pp = 0.75;
    const std::vector<double> probabilities = boost::assign::list_of(pp)(1-pp);
    std::ofstream fout("mn_2D_parameters.in");
    fout << "N = " << NN << std::endl;
    fout << "p = " << pp << std::endl;

    /* create 2D multinomial_distribution object */
    mn_type mn_rnd(NN, probabilities);

    /* number of samples to draw from distribution */
    const int number_trials = 1e3;

    /* sample space */
    std::vector<int_type> sample_space(NN + 1, -1.0);
    for (int_type ii = 0; ii < sample_space.size(); ii++) {
        sample_space.at(ii) = ii;
    }

    /* create CDF object */
    probability::CDF_discrete mn_cdf(sample_space, "mn_2D_cdf.dat");

    /* repeatedly add samples to the CDF */
    for (int ii = 0; ii < number_trials; ii++) 
        mn_cdf.update(mn_rnd(urng).at(0));

    std::cout << "... CDF generated" << std::endl;
    
}


/*************************************************************************/


void unit_test::unit_test__sum_random_vector_elements() {

    std::cout << "testing sum of elements in random vector..." << std::endl;

    /* PRNG engine */
    boost::random::mt19937 urng;

    /* parameters of multinomial_distribution */
    const int_type NN = 10;
    const std::vector<double> probabilities = boost::assign::list_of(0.2)(0.2)(0.2)(0.4);

    /* create multinomial_distribution object */
    mn_type mn_rnd(NN, probabilities);
    std::cout << "generated " << probabilities.size() << "D multinomial distribution" << std::endl;

    /* number of samples to draw from distribution */
    const int number_trials = 1e2;

    /* sample multinomial distribution */
    for (int ii = 0; ii < number_trials; ii++) {

        std::vector<int_type> random_vector = mn_rnd(urng);

        int_type sum = static_cast<int_type>(0);
        for (int ii = 0; ii < random_vector.size(); ii++) {
            sum += random_vector.at(ii);
        }

        assert(sum == NN);

    }

    std::cout << "... finished testing sum of elements in random vector" << std::endl;

}

/*************************************************************************/


void unit_test::unit_test__3Dmultinomial() {

    std::cout << "sampling from 3D multinomial ... " << std::endl;

    /* PRNG engine */
    boost::random::mt19937 urng;

    /* parameters of multinomial_distribution */
    const int_type NN = 20;
    const std::vector<double> probabilities = boost::assign::list_of(0.8)(0.1)(0.1);
    std::ofstream fin("mn_3D_parameters.in");
    fin << "N = " << NN << std::endl;
    for (int ii = 0; ii < probabilities.size(); ii++)
        fin << "p" << ii << " = " << probabilities.at(ii) << std::endl;

    /* create multinomial_distribution object */
    mn_type mn_rnd(NN, probabilities);
    std::cout << "generated " << probabilities.size() << "D multinomial distribution" << std::endl;

    /* number of samples to draw from distribution */
    const int number_trials = 1e5;

    /* open file to write samples to */
    std::ofstream fout("mn_3D_samples.dat");
    
    /* sample multinomial distribution */
    for (int ii = 0; ii < number_trials; ii++) {

        std::vector<int_type> random_vector = mn_rnd(urng);

        for (int ii = 0; ii < random_vector.size(); ii++) {
            fout << random_vector.at(ii) << " ";
        }
        fout << std::endl;

    }

    std::cout << "... finished sampling from 3D multinomial" << std::endl;

}

/*************************************************************************/


void unit_test::unit_test__4Dmultinomial() {

    std::cout << "sampling from 4D multinomial ... " << std::endl;

    /* PRNG engine */
    boost::random::mt19937 urng;

    /* parameters of multinomial_distribution */
    const int_type NN = 10;
    const std::vector<double> probabilities = boost::assign::list_of(0.6)(0.1)(0.1)(0.2);
    std::ofstream fin("mn_4D_parameters.in");
    fin << "N = " << NN << std::endl;
    for (int ii = 0; ii < probabilities.size(); ii++)
        fin << "p" << ii << " = " << probabilities.at(ii) << std::endl;

    /* create multinomial_distribution object */
    mn_type mn_rnd(NN, probabilities);
    std::cout << "generated " << probabilities.size() << "D multinomial distribution" << std::endl;

    /* number of samples to draw from distribution */
    const int number_trials = 1e5;

    /* open file to write samples to */
    std::ofstream fout("mn_4D_samples.dat");

    /* sample multinomial distribution */
    for (int ii = 0; ii < number_trials; ii++) {

        std::vector<int_type> random_vector = mn_rnd(urng);

        for (int ii = 0; ii < random_vector.size(); ii++) {
            fout << random_vector.at(ii) << " ";
        }
        fout << std::endl;

    }

    std::cout << "... finished sampling from 4D multinomial" << std::endl;

}


/*************************************************************************/


void unit_test::unit_test__3Dmultinomial_edgeCase() {

    std::cout << "sampling from 3D multinomial (edge case)... " << std::endl;

    /* PRNG engine */
    boost::random::mt19937 urng;

    /* parameters of multinomial_distribution */
    const int_type NN = 5;
    const std::vector<double> probabilities = boost::assign::list_of(0.0)(0.0)(1.0);
    std::ofstream fin("mn_3D_edgeCase_parameters.in");
    fin << "N = " << NN << std::endl;
    for (int ii = 0; ii < probabilities.size(); ii++)
        fin << "p" << ii << " = " << probabilities.at(ii) << std::endl;

    /* create multinomial_distribution object */
    mn_type mn_rnd(NN, probabilities);
    std::cout << "generated " << probabilities.size() << "D multinomial distribution" << std::endl;

    /* number of samples to draw from distribution */
    const int number_trials = 1e5;

    /* sample multinomial distribution */
    for (int ii = 0; ii < number_trials; ii++) {

        std::vector<int_type> random_vector = mn_rnd(urng);

        assert(random_vector.at(0) == 0);
        assert(random_vector.at(1) == 0);
        assert(random_vector.at(2) == NN);

//        for (int ii = 0; ii < random_vector.size(); ii++) {
//            std::cout << random_vector.at(ii) << " ";
//        }
//        std::cout << std::endl;

    }

    std::cout << "... finished sampling from 3D multinomial (edge case)" << std::endl;

}



/*************************************************************************/

void unit_test::unit_test__3Dmultinomial_zeroN() {

    std::cout << "sampling from 3D multinomial (N = 0)... " << std::endl;

    /* PRNG engine */
    boost::random::mt19937 urng;

    /* parameters of multinomial_distribution */
    const int_type NN = 0;
    const std::vector<double> probabilities = boost::assign::list_of(0.1)(0.1)(0.8);
    std::ofstream fin("mn_3D_zeroN_parameters.in");
    fin << "N = " << NN << std::endl;
    for (int ii = 0; ii < probabilities.size(); ii++)
        fin << "p" << ii << " = " << probabilities.at(ii) << std::endl;

    /* create multinomial_distribution object */
    mn_type mn_rnd(NN, probabilities);
    std::cout << "generated " << probabilities.size() << "D multinomial distribution" << std::endl;

    /* number of samples to draw from distribution */
    const int number_trials = 1e5;

    /* sample multinomial distribution */
    for (int ii = 0; ii < number_trials; ii++) {

        std::vector<int_type> random_vector = mn_rnd(urng);

        assert(random_vector.at(0) == 0);
        assert(random_vector.at(1) == 0);
        assert(random_vector.at(2) == 0);

//        for (int ii = 0; ii < random_vector.size(); ii++) {
//            std::cout << random_vector.at(ii) << " ";
//        }
//        std::cout << std::endl;

    }

    std::cout << "... finished sampling from 3D multinomial (N = 0)" << std::endl;

}



