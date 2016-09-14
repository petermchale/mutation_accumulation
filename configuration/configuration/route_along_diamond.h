#ifndef ROUTE_ALONG_DIAMOND_H
#define	ROUTE_ALONG_DIAMOND_H

/*************************************************************************/

namespace monte_carlo {

    /**
     * indicates which single-mutant gave rise to first K-fold mutant in each sub-population\n
     * assumes discrete branching process with K=2 (Branching_Discrete_Diamond) \n
     * in a discrete branching process, both single-mutant species may generate a double-mutant stem cell in the same cell cycle\n
     * generalize to handle arbitrary number of K-1-fold mutants \n
     */
    template <class population_type>
    class Diamond_Routes {
    private:

        /** 
         * _a_yielded_ab(pop) indicates whether first ab arose from mutation of an "a" in sub-population pop \n
         * uses deque to handle bool template parameter \n
         * see: http://stackoverflow.com/questions/3458856/overloaded-operator-on-template-class-in-c-with-const-nonconst-versions
         */
        std::deque<bool> _a_yielded_ab;
        std::deque<bool> _b_yielded_ab;
        std::deque<bool> _first_ab_generated;

    public:

        /**
         * custom constructor
         */
        explicit Diamond_Routes(const int &number_sub_pops_) {

            _a_yielded_ab = std::deque<bool> (number_sub_pops_, false);
            _b_yielded_ab = std::deque<bool> (number_sub_pops_, false);
            _first_ab_generated = std::deque<bool> (number_sub_pops_, false);

        }

        /**
         * update routes
         */
        void update(
                const int &pop,
                const population_type &number_ab_generated_by_a,
                const population_type &number_ab_generated_by_b) {

            if (!_first_ab_generated.at(pop)) { // if ab has not yet arisen ...
                if (number_ab_generated_by_a > static_cast<population_type> (0)) // if one or more "a" stem cells have generated ab(s) ...
                {
                    _a_yielded_ab.at(pop) = true; // indicate that "a" has generated "ab"
                    _first_ab_generated.at(pop) = true; // indicate that "ab" has been generated
                }
                if (number_ab_generated_by_b > static_cast<population_type> (0)) // if one or more "b" stem cells have generated ab(s) ...
                {
                    _b_yielded_ab.at(pop) = true; // indicate that "b" has generated "ab"
                    _first_ab_generated.at(pop) = true; // indicate that "ab" has been generated
                }
            }

        }

        /**
         * query route a
         */
        const bool a_yielded_ab(const Pop &pop) const {

            return _a_yielded_ab.at(pop.value());
        }

        /**
         * query route b
         */
        const bool b_yielded_ab(const Pop &pop) const {

            return _b_yielded_ab.at(pop.value());
        }


    };


}

#endif	/* ROUTE_ALONG_DIAMOND_H */

