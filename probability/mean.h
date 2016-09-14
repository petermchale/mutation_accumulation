#ifndef MEAN_CLASS_H
#define	MEAN_CLASS_H

#include <cassert> // assert

#include "moment.h" // probability::Moment 

/*************************************************************************/

namespace probability {

    /**
     * calculate mean from ensemble of samples \n
     * may decorate this class to estimate convergence (Joshi's book pp 72 - 81)\n
     * use template arguments instead of Joshi's wrapper class to make decorator\n
     */
    template <class sample_type>
    class Mean : public Moment<sample_type> {
    private:

        sample_type _running_sum;
        long long int _number_trials;
        
    public:

        /** 
         * custom ctor 
         */
        explicit Mean() : _running_sum(static_cast<sample_type>(0)), _number_trials(0) {
            
        }
        
        /**
         * update estimate of moment
         */
        virtual void update(const sample_type &sample) {

            _running_sum += sample; 
            _number_trials++; 

        }

        /**
         * get moment
         */
        virtual const double value() const {
            
            if (_number_trials > 0) 
                return _running_sum/static_cast<double>(_number_trials);
            else 
                assert(false);
        }

        /**
         * return number trials
         */
        virtual const long long int trials() const { 
            
            return _number_trials;
        }

    };






}

#endif	/* MEAN_CLASS_H */


