# mutation_accumulation
C++ codebase to simulate the stochastic accumulation of mutations in individual stem cells. This code was used to produce the results reported in [this cancer systems biology paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003802)

Templates are a feature of the C++ programming language that allows classes to operate with generic types. This allows a function or class to work on many different data types without being rewritten for each one. 

Example
============
* download the repository and make sure it is called 'mutation_accumulation'
* navigate to the example sub-directory
* make your C++ compiler aware of the path of the mutation_accumulation directory 
* compile main_branching_CDF_trajs.cpp
* run the executable in the example directory (you'll find pre-existing output in the data directory)
* run plot_trajs.m in Matlab to see a time course tracking the number of wild-type, single-mutant, and double-mutant stem cell in a stochastic simulation prior to the appearance of the first triple-mutant stem cell
* run plot_cdf.m in Matlab to see the cumulative probability that the first triple-mutant stem cell has arisen as a function of time (measured in units of the average time it takes a stem cell to divide)
* look at this [Jupyter Notebook]() to see a mathematical analysis of the Monte Carlo simulation


