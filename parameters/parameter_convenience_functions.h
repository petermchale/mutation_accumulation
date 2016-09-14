#ifndef PARAMETER_CONVENIENCE_FUNCTIONS_H
#define	PARAMETER_CONVENIENCE_FUNCTIONS_H

#include <sstream> // std::stringstream

namespace parameters {

    /**
     * convert a model parameter value into a string\n
     *
     * prefer non-member function to member function\n
     * eg. Symmetry "is-a" Parameter_base<double> (Item 32)
     * 
     */
    template <class value_type>
    const std::string get_str(const Parameter_base<value_type> & parameter) {

        std::stringstream ss; // create a stringstream
        ss << parameter.value(); // add value to the stream
        return ss.str(); // return a string with the contents of the stream
    }


} // end namespace parameters

#endif	/* PARAMETER_CONVENIENCE_FUNCTIONS_H */

