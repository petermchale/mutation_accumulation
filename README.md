# mutation_accumulation
C++ template programming framework to perform Monte Carlo simulations of the accumulation of mutations in individual stem cells. 

###Scientific use
You'll find a brief mathematical analysis in this [Jupyter Notebook]() of the Monte Carlo simulation in the `example` directory. `mutation_accumulation` was written primarily to produce the results reported in [this cancer systems biology paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003802). 


###Getting started
* download the repository and rename it to `mutation_accumulation`
* add the directory in which `mutation_accumulation` resides to the list of include directories that your C++ compiler searches
* download the Boost library and similarly make your compiler aware of its location

```unix
g++ -I<path to boost> -I<path to mutation_accumulation> main_branching_CDF_trajs.cpp
```

* navigate to the example sub-directory
* compile main_branching_CDF_trajs.cpp
* run the executable in the example directory (you'll find pre-existing output in the data directory)
* run plot_trajs.m in Matlab to see a time course tracking the number of wild-type, single-mutant, and double-mutant stem cell in a stochastic simulation prior to the appearance of the first triple-mutant stem cell
* run plot_cdf.m in Matlab to see the cumulative probability that the first triple-mutant stem cell has arisen as a function of time (measured in units of the average time it takes a stem cell to divide)

###Template programming
Templates are a feature of the C++ programming language that allows classes to operate with generic types. This allows a function or class to work on many different data types without being rewritten for each one. 
```C++
    typedef long long int population_type;
    typedef monte_carlo::Branching_Discrete<population_type> Configuration_Policy;

    typedef Configuration_Policy::time_t time_type;
    typedef probability::Notify_NonNegative_BoundedAbove<time_type> Notification_Policy;
    typedef probability::CDF<Notification_Policy> Histogram_Policy;

    monte_carlo::Calculate_Histogram_Trajs<Histogram_Policy, Configuration_Policy, monte_carlo::Raw_Data_Null,  monte_carlo::Read_NonHomeostasis_Policy>::implement();
```

