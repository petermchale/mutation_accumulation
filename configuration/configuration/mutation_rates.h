#ifndef MUTATION_RATES_H
#define	MUTATION_RATES_H

#include <mutation_accumulation/array/array1D.h> // array::Array1D

/*************************************************************************/

namespace monte_carlo {

    namespace MutationRates_namespace {

        typedef double data_type;

        /**
         * 1D array to hold mutation rates \n
         * making Array1D a private data member (instead of inheriting it) would be better Fowler p352
         */
        class MutationRates : public array::Array1D<data_type> {
        public:

            typedef data_type data_t;

        public:

            /**
             * custom constructor
             */
            explicit MutationRates(const std::vector<data_type> &uu)
            : array::Array1D<data_type>(uu.size()) {

                data = uu;
            }

            /**
             * must declare default constructor: \n
             * compiler provides default constructor only if no custom constructors exist.
             */
            explicit MutationRates() : array::Array1D<data_type>() {

            }

            /**
             * custom constructor
             */
            explicit MutationRates(const int &size_)
            : array::Array1D<data_type>(size_) {

            }

            /**
             * number of mutation rates
             */
            const int size() const {

                return (this->dim);
            }

        };

    }

    using MutationRates_namespace::MutationRates;



}

#endif	/* MUTATION_RATES_H */

