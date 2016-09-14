#ifndef SAMPLE_SPACE_H
#define	SAMPLE_SPACE_H

#include <mutation_accumulation/utility/grid.h> // grid::make_uniform_grid

/*************************************************************************/

namespace probability {

    /**
     * create bernoulli sample space = (false, true) 
     */
    const std::vector<bool> make_bernoulli_sample_space() {

        std::vector<bool> sample_space(2);
        sample_space.at(0) = false;
        sample_space.at(1) = true;

        return sample_space;
    }

    /**
     * create uniform sample space of arbitrary type starting at zero\n
     * deprecated in favor of Make_Uniform_Sample_Space\n
     * retained for backward compatibility\n
     */
    template <class sample_type>
    const std::vector<sample_type> make_uniform_sample_space(
    const int &suggested_number_sample_points,
    const sample_type &span_of_sample_space) {

        std::cerr << "make_uniform_sample_space is deprecated in favor of Make_Uniform_Sample_Space" << std::endl;

        return grid::make_uniform_grid(suggested_number_sample_points, static_cast<sample_type>(0), span_of_sample_space);
    }

    /** 
     * create uniform sample space of arbitrary type starting at zero\n
     */
    template <class sample_type>
    class Make_Uniform_Sample_Space {
    public:

        typedef sample_type sample_t;

        static const std::vector<sample_type> build(
                const int &suggested_number_sample_points,
                const sample_type &sample_space_lower,
                const sample_type &sample_space_upper) {

            return grid::make_uniform_grid(suggested_number_sample_points, sample_space_lower, sample_space_upper);

        }

    };

    /** 
     * create logarithmic sample space of arbitrary type\n
     */
    template <class sample_type>
    class Make_Logarithmic_Sample_Space {
    public:

        typedef sample_type sample_t;

        static const std::vector<sample_type> build(
                const int &suggested_number_sample_points,
                const sample_type &sample_space_lower,
                const sample_type &sample_space_upper) {

            return grid::make_logarithmic_grid(suggested_number_sample_points, sample_space_lower, sample_space_upper);

        }

    };

}

#endif	/* SAMPLE_SPACE_H */


