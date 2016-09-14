#ifndef PARAMETER_BASE_H
#define	PARAMETER_BASE_H

/**
 * this namespace includes wrapper classes that hold values of model parameters (strong typing) \n
 * such classes make application code easier to read\n
 * such classes guard against passing function/method arguments in the wrong order\n
 * see Scott Meyer's book for rationale
 */
namespace parameters {

    /**
     * Abstract base class for model parameters \n
     * Should make protected data members private (Item 22)
     */
    template <class value_type>
    class Parameter_base {
        
    protected:

        value_type value_;
        value_type value_min_;
        value_type value_max_;

    protected:

        /**
         * ctor
         */
        explicit Parameter_base(
                const value_type &value__,
                const value_type &value_min__,
                const value_type &value_max__) :
        value_(value__), value_min_(value_min__), value_max_(value_max__) {

        };

    public:
        
        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers
         */
        virtual ~Parameter_base() {

        }


        /**
         * get value of model parameter
         */
        const value_type value() const {

            return value_;
        }

        
    };

} // end namespace parameters

#endif	/* PARAMETER_BASE_H */

