#ifndef OBSERVER_H
#define	OBSERVER_H

#include <list> // std::list

#include <boost/shared_ptr.hpp> // boost::shared_ptr
#include <boost/utility.hpp> // boost::noncopyable

/*************************************************************************/

/**
 * design patterns
 */
namespace patterns {

    /* note that boost::signals can do Observer pattern, but boost::signals library must be linked as opposed to simply include'd */

    
    /**
     * see Item 6 for use of boost::noncopyable\n
     * private inheritance means is-implemented-in-terms-of (Item 39)\n
     * C++ guarantees that the ctor and dtor of Observer and Subject call the ctor and dtor of boost::noncopyable (Item 30)
     */


    /**
     * abstract base class that defines the Observer interface\n
     * cannot instantiate this class because of pure virtual member function\n
     * this class is uncopyable to prevent the shallow copying of member references/pointers (in derived classes)
     */
    class Observer : private boost::noncopyable {
        
    public:

        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers\n
         * 
         */
        virtual ~Observer() {

        }

        /**
         * update Observer\n
         */
        virtual void update() = 0; // re-declare and implement in derived classes

    };

    

    /**
     * base class that defines the Subject interface\n
     * cannot instantiate this class because ctor is protected\n
     * this class is uncopyable to prevent shallow copying of member references/pointers
     *
     */
    class Subject : private boost::noncopyable {
    private:

        /**
         * things in a container cannot be references\n
         * thus cannot use "Observer &"
         *
         * things in a container must be assignable (i.e. not const)\n
         * thus cannot use "Observer * const", even though the pointers are not going to be changed
         */
        typedef std::list<Observer * > list_t;
        list_t _observers;

    protected:

        /**
         * default ctor
         */
        Subject() {

            /* _observers' default ctor, std::list(), should make an empty list with no content and a size of zero */
            if (!_observers.empty()) {
                
                std::cerr << "observer list is not empty!" << std::endl;
                assert(false);

            }
        }


        /**
         * notify observers\n
         */
        void notify() const {

            typedef list_t::const_iterator list_it_t;
            for (list_it_t it = _observers.begin(); it != _observers.end(); it++) {
                (*it)->update();
            }
        }


    public:
        
        /**
         * virtual destructor\n
         * ensures that objects can be deleted properly using interface base class pointers
         */
        virtual ~Subject() {

        }

        /**
         * attach Observer\n
         */
        void attach(Observer * const observer_ptr) {

            _observers.push_back(observer_ptr);
            
        }

        
        /**
         * detach Observer\n
         * prevents dangling references
         */
        void detach(Observer * const observer_ptr) {

            _observers.remove(observer_ptr); // remove element if pointers compare equal
            
        }


    };


    /**
     * factory function that dynamically allocates an Observer AND attaches it to a Subject\n
     * Items 13, 18 of Meyer's book\n
     * factory assumes that Observer writes a log file every observer_divisor notifications\n
     * could make more generic by lumping all observer-dependent arguments into a class
     */
    template <class Obs_t, class Sub_t>
    boost::shared_ptr<Obs_t> createAttachedObserver(Sub_t &subject, const std::string &logFileName, const int &observer_divisor) {

        boost::shared_ptr<Obs_t> observer_ptr(new Obs_t(subject, logFileName, observer_divisor));

        observer_ptr->attach_to_subject();

        return observer_ptr;
    }


    /**
     * factory function that dynamically allocates a Subject\n
     * Items 13, 18 of Meyer's book\n
     * factory assumes the ctor specified by CDF \n
     */
    template <class Sub_t, class sample_type>
    boost::shared_ptr<Sub_t> createSubject(const std::vector<sample_type> &sample_space, const std::string &fileName) {

        boost::shared_ptr<Sub_t> subject_ptr(new Sub_t(sample_space, fileName));

        return subject_ptr;
    }


} 

#endif	/* OBSERVER_H */


