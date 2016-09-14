#ifndef FITNESS_H
#define	FITNESS_H

/*************************************************************************/

namespace monte_carlo {

    namespace Fitness_namespace {

        typedef double data_type;

        /**
         * 1D array to hold fitness for each stage (except the last) \n
         */
        class Fitness {
        public:

            typedef data_type data_t;

        private:

            std::vector<data_type> _ww;

        public:

            /**
             * custom constructor
             */
            explicit Fitness(const std::vector<data_type> &ww) : _ww(ww) {

            }

            /**
             * must declare default constructor: \n
             * compiler provides default constructor only if no custom constructors exist.
             */
            explicit Fitness() {

                _ww = std::vector<data_type > ();
            }

            /**
             * custom constructor
             */
            explicit Fitness(const int &size_) {

                _ww = std::vector<data_type > (size_);
            }

            /**
             * number of fitness values
             */
            const int size() const {

                return _ww.size();
            }

            /**
             * read particular fitness value
             */
            const data_type at(const int &ii) const {

                return _ww.at(ii);
            }

        };

    }

    using Fitness_namespace::Fitness;



}

#endif	/* FITNESS_H */

