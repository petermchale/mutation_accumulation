#ifndef CREATE_CONFIGURATION_H
#define	CREATE_CONFIGURATION_H

#include <mutation_accumulation/utility/configuration_traits.h> // configuration_types::branching_type etc
#include <mutation_accumulation/simulation/read_policy.h> // Read_Homeostasis_Policy<Configuration_type>

/*************************************************************************/

namespace monte_carlo {

    /** 
     * create instance of a moran configuration (renewal bias must equal 1/2) 
     */
    template <class Configuration_type>
    const Configuration_type create_configuration(
            const configuration_categories::moran_category,
            const Read_Homeostasis_Policy<Configuration_type> &read_homeostasis_policy) {

        return Configuration_type(
                read_homeostasis_policy.get_population(),
                read_homeostasis_policy.getUu(), 
                read_homeostasis_policy.getSymmetry(), 
                read_homeostasis_policy.getTime_grid());

    }

    /** 
     * create instance of a branching configuration, with arbitrary renewal bias
     */
    template <class Configuration_type>
    const Configuration_type create_configuration(
            const configuration_categories::branching_category,
            const Read_Policy_Base<Configuration_type> &read_policy) {

        return Configuration_type(
                read_policy.get_population(),
                read_policy.getUu(), 
                read_policy.getSymmetry(), 
                read_policy.get_SymmetricRenewal(), 
                read_policy.getTime_grid());

    }

}

#endif	/* CREATE_CONFIGURATION_H */

