/* 
 * File:   multinomial_distribution.h
 * Author: petermchale
 *
 * Created on August 14, 2012, 2:18 PM
 */

#ifndef MULTINOMIAL_DISTRIBUTION_H
#define	MULTINOMIAL_DISTRIBUTION_H

#include <vector> //std::vector

#include <boost/random/binomial_distribution.hpp> // boost::random::binomial_distribution

/*************************************************************************/

namespace mutation_accumulation {

    namespace random {

        /**
         * sample from the multinomial distribution using the "conditional method"; see \n
         * Devroye's book (1993) p558;\n
         * "The computer generation of multinomial random variates, CS Davis, Comp Stat & Data Analysis, 16 (1993), p205-217\n
         * \n
         *
         * For multiple template type parameters, all parameters after the first default argument must have default parameters.\n
         * When declaring a template class object with default type parameters, omit the parameters to accept the default argument.\n
         * If there are no nondefault parameters, do not omit the empty angle brackets.
         */
        template<class IntType = int, class ProbType = double>
        class multinomial_distribution {
        private:

            IntType _N;
            std::vector<ProbType> _probabilities;

        public:

            /**
             * Construct a @c multinomial_distribution object.\n
             * @c N and @c probabilities are the parameters of the distribution.\n
             * Requires: N >=0 && \sum_i probabilities_i = 1
             */
            explicit multinomial_distribution(const IntType& N_arg, const std::vector<ProbType>& probabilities_arg)
            : _N(N_arg), _probabilities(probabilities_arg) {

            }

            /**
             * Returns a random vector distributed according to the multinomial distribution. \n
             * The elements of the random vector sum to N.
             */
            template<class URNG>
            const std::vector<IntType> operator()(URNG& urng) const {

                if (_N == 0) {

                    std::vector<IntType> random_vector(_probabilities.size(), static_cast<IntType> (0));

                    return random_vector;

                } else if (_N > 0) {

                    IntType N_bar = _N;
                    ProbType qq = static_cast<ProbType> (1);
                    std::vector<IntType> random_vector(_probabilities.size(), static_cast<IntType> (-1));

                    for (int index = 0; index < _probabilities.size() - 1; index++) {

                        typedef boost::random::binomial_distribution<IntType, ProbType> binomial_distribution_type;
                        const ProbType p_bar = _probabilities.at(index) / qq;
                        random_vector.at(index) = binomial_distribution_type(N_bar, p_bar)(urng);

                        N_bar -= random_vector.at(index);
                        qq -= _probabilities.at(index);

                    }

                    random_vector.back() = N_bar;

                    return random_vector;

                } else {

                    std::cerr << "multinomial_distribution: _N = " << _N << std::endl;
                    exit(1);

                }

            }

            /**
             * vector of probabilities that defines multinomial distribution
             */
            const std::vector<ProbType> get_probabilities() const {

                return _probabilities;
            }

            /**
             * N parameter of multinomial distribution
             */
            const IntType get_N() const {

                return _N;
            }



        };
    }
}

#endif	/* MULTINOMIAL_DISTRIBUTION_H */

