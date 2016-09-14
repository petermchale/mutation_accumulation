#ifndef MOMENT_BASE_H
#define	MOMENT_BASE_H

/*************************************************************************/

namespace probability {

    /**
     * abstract base class that defines the interface that moment sub-classes must conform to\n
     * a moment is one of: mean, variance, skewness, kurtosis, etc
     */
    template <class sample_type>
    class Moment {
    public:

        typedef sample_type sample_t;
        typedef double Results_t;
        
    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * this function is public so that one can call delete through a base-class pointer
         */
        virtual ~Moment() {

        }

        /**
         * update estimate of moment
         */
        virtual void update(const sample_type &sample) = 0;

        /**
         * get moment
         */
        virtual const Results_t value() const = 0;

        /**
         * return number trials
         */
        virtual const long long int trials() const = 0;

    };






} 

#endif	/* MOMENT_BASE_H */


