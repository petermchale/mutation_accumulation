#ifndef FILES_H
#define	FILES_H

#include <boost/shared_ptr.hpp> // boost::shared_ptr
#include <fstream> // std::ifstream

/*************************************************************************/

namespace monte_carlo {

    /**
     * write an empty file called "done"\n
     * use "find . -name done" to find sub-directories in which simulation has completed
     */
    inline void done() {

        std::string fileName = "done";
        std::fstream fout(fileName.c_str(), std::ios_base::out); // file closed when std::fstream object goes out of scope (RAII)

        if (!fout.is_open()) {
            std::cerr << "cannot open " << fileName << std::endl;
            assert(false);
        }

    }

    /** 
     * opens file for input and checks that file opened correctly 
     */
    inline const boost::shared_ptr<std::ifstream> open_file_for_input(const std::string &fileName) {

        // file will be closed when last shared_ptr to following dynamically allocated ifstream object goes out of scope (RAII)
        boost::shared_ptr<std::ifstream> ifstream_ptr(new std::ifstream(fileName.c_str()));

        if (!ifstream_ptr->is_open()) {
            std::cerr << "cannot open " << fileName << std::endl;
            assert(false);
        }

        return ifstream_ptr;

    }

    /** 
     * opens file for output and checks that file opened correctly 
     */
    inline const boost::shared_ptr<std::ofstream> open_file_for_output(const std::string &fileName) {

        // file will be closed when last shared_ptr to following dynamically allocated ofstream object goes out of scope (RAII)
        boost::shared_ptr<std::ofstream> ofstream_ptr(new std::ofstream(fileName.c_str()));

        if (!ofstream_ptr->is_open()) {
            std::cerr << "cannot open " << fileName << std::endl;
            assert(false);
        }

        return ofstream_ptr;

    }



}



#endif	/* FILES_H */

