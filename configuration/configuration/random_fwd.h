#ifndef RANDOM_FWD_H
#define	RANDOM_FWD_H

#include <boost/random/mersenne_twister.hpp> // boost::random::mt19937

/*************************************************************************/

namespace monte_carlo {

    typedef boost::random::mt19937 base_generator_type; 

}

namespace patterns {

    typedef boost::random::mt19937 base_generator_type; 

}

#endif	/* RANDOM_FWD_H */

