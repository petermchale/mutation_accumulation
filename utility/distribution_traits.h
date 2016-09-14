#ifndef DISTRIBUTION_TRAITS_H
#define	DISTRIBUTION_TRAITS_H

/*************************************************************************/

namespace distribution_types {

    /*
     * These are empty types, used to distinguish different types of probability distribution.\n
     * User-defined probability distribution types must publicly typedef "category" to alias one of these empty types
     */

    /**
     *  marks probability mass functions with boolean sample space
     */
    struct pmf_bool_type {
    };

    /**
     *  marks cumulative distribution functions (discrete or continuous sample space)
     */
    struct cdf_type {
    };

}



#endif	/* DISTRIBUTION_TRAITS_H */


