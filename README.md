# BioCCP
*The Coupon Collector's Problem for Combinatorial Biotechnology*

During the **combinatorial engineering of biosystems** such as proteins, genetic circuits and genomes, diverse libraries are generated by assembling and recombining modules. the variants with the optimal phenotype are selected with screening techniques. However, when the number of available modules to compose a biological designs increases, a combinatorial explosion of design possibilities arises, allowing only for a part of the libary to be analyzed. In this case, it is important for a researcher to get insight in which (minimal) sample size sufficiently covers the design space, i.e. what is the **minimal number of designs so that all modules are observed at least once**.

BioCCP contains functions for calculating minimum sample sizes and related statistics:


Function name    | Short description
---------------- | -----------------
`expectation_minsamplesize`        | Calculated expectation of required number of designs to observe all modules at least *m* times
`std_minsamplesize`      | Calculated standard deviation on required number of designs to observe all modules at least *m* times
`success_probability`         | Calculates the chance that a given number of designs is less than or equal to the required number of designs to sea each module at least *m* times
`expectation_fraction_collected` | Returns for a given sample size the expected fraction of modules in the design space observed
`prob_occurence_module` | Calculates for a module with specified module probability *p*, the chance that the module has occured *j* times when a given number of designs has been collected
 

For more info about the implementation of the functions, please consult the docs.

The Pluto notebook **BioCCP Interactive Notebook.jl** provides an interactive illustration of these functionalities and allows to make calculations for specific settings. 

## Getting started

  

