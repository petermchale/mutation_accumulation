
#include <mutation_accumulation/simulation/calculate_histogram_trajs.h> 
#include <mutation_accumulation/configuration/configuration/branching_discrete.h> 
#include <mutation_accumulation/probability/notification_policy.h> 
#include <mutation_accumulation/probability/cdf.h> 
#include <mutation_accumulation/simulation/raw_data.h>

/*************************************************************************/

int main() {

    typedef long long int population_type;
    typedef monte_carlo::Branching_Discrete<population_type> Configuration_Policy;

    typedef Configuration_Policy::time_t time_type;
    typedef probability::Notify_NonNegative_BoundedAbove<time_type> Notification_Policy;
    typedef probability::CDF<Notification_Policy> Histogram_Policy;

    monte_carlo::Calculate_Histogram_Trajs<Histogram_Policy, Configuration_Policy, monte_carlo::Raw_Data_Null, monte_carlo::Read_NonHomeostasis_Policy>::implement();
}



