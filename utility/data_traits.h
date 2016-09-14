#ifndef DATA_TRAITS_H
#define	DATA_TRAITS_H

/*************************************************************************/

namespace data_types {

    /*
     *  These are empty types, used to distinguish different types of data.  The
     *  distinction is not made by what they contain, but simply by what they
     *  are.  Different underlying algorithms can then be used based on the
     *  different operations supported by different data types.
     */

    /**
     *  marks discrete built-in data types 
     */
    struct discrete_type {
    };

    /**
     *  marks continuous built-in data types 
     */
    struct continuous_type {
    };

    /**
     *  This class does nothing but define nested typedefs.  The general
     *  version imposes the requirement that any user-defined data type must contain a nested typedef
     *  named "category" that identifies whether the data type is discrete or continuous.  
     */
    template <class data_type>
    struct data_traits {
        typedef typename data_type::category category;
    };

    /* 
     *  Specialized versions of traits class 
     *  indicates whether commonly used built-in data types are discrete or continuous.
     */
    /** 
     * int is discrete 
     */
    template <>
    struct data_traits<int> {
        typedef discrete_type category;
    };

    /** 
     * long long int is discrete 
     */
    template <>
    struct data_traits<long long int> {
        typedef discrete_type category;
    };
    /** 
     * bool is discrete 
     */
    template <>
    struct data_traits<bool> {
        typedef discrete_type category;
    };

    /** 
     * double is continuous
     */
    template <>
    struct data_traits<double> {
        typedef continuous_type category;
    };
}



#endif	/* DATA_TRAITS_H */


