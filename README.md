# BioCCP
*The Coupon Collector's Problem for Combinatorial Biotechnology*

## Intro
During the **combinatorial engineering of biosystems** such as proteins, genetic circuits and genomes, diverse libraries are generated by assembling and recombining modules. The variants with the optimal phenotype are selected with screening techniques. However, when the number of available modules to compose a biological designs increases, a combinatorial explosion of design possibilities arises, allowing only for a part of the libary to be analyzed. In this case, it is important for a researcher to get insight in which (minimal) sample size sufficiently covers the design space, i.e. what is the **minimal number of designs so that all modules are observed at least once**.


## Functions
BioCCP contains functions for calculating minimum sample sizes and related statistics:

Function name    | Short description
---------------- | -----------------
`expectation_minsamplesize`        | Calculated expectation of required number of designs to observe all modules at least *m* times
`std_minsamplesize`      | Calculated standard deviation on required number of designs to observe all modules at least *m* times
`success_probability`         | Calculates the chance that a given number of designs is less than or equal to the required number of designs to sea each module at least *m* times
`expectation_fraction_collected` | Returns for a given sample size the expected fraction of modules in the design space observed
`prob_occurence_module` | Calculates for a module with specified module probability *p*, the chance that the module has occured *j* times when a given number of designs has been collected
 

For more info about the implementation of the functions, please consult the docs.

## Pluto notebook

The Pluto notebook provides an interactive illustration of all functions in BioCCP and assembles a report for your specific design set-up. 

Inputs for generating the report:
- total number of modules in the design space, 
- the number of modules per design, 
- the number of complete sets of modules to collect, and
- probability distribution of the modules (custom or by assuming a specific distribution)
   - *When exact probabilities are known:* A custom module probability/abundance vector 
   - *When probabilities unknown:
      - A max        

Using this information, a report for sample size determination is created using the [functions](## Functions) described above. The report contains the following sections:
- Module probabilities
- Minimum sample size required
- Success probability
- Expected saturation
- Occurence of a specific module

* *the notebook allows to enter a custom module probabilty/abundance vector or to assume that the module probabilities obey one of the following distributions:
- 


## Getting started

#### Launch Pluto notebook from Browser 

Launch the Pluto notebook directly from your browser using Binder (no installation of Julia/packages required): [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kirstvh/PlutoNotebooks/main?urlpath=pluto/open?path=/home/jovyan/notebooks/BioCCP_Interactive_Notebook.jl)

#### Execute functions in Julia

(1) [Install Julia](https://julialang.org/downloads/) 

(2) Install BioCCP in the Julia REPL (hit ] to enter package mode):

    ] add https://github.com/kirstvh/BioCCP

(3) For using the [Pluto notebook](BioCCP/notebook/BioCCP_Interactive_Notebook.jl):

 In the Julia REPL, hit the following command to install the additional packages [Pluto](https://github.com/fonsp/Pluto.jl), [PlutoUI](https://github.com/fonsp/PlutoUI.jl) and [Plots](https://github.com/JuliaPlots/Plots.jl) 
  
    ] add Pluto, PlutoUI, Plots

 Then start Pluto in the Julia REPL:

    using Pluto; Pluto.run()

## References
Implementation of formula's was based on:

> Doumas, A. V., & Papanicolaou, V. G. (2016). *The coupon collector’s problem revisited: generalizing the double Dixie cup problem of Newman and Shepp.* ESAIM: Probability and Statistics, 20, 367-399. doi: 	https://doi.org/10.1051/ps/2016016

> Boneh, A., & Hofri, M. (1997). *The coupon-collector problem revisited—a survey of engineering problems and computational methods.* Stochastic Models, 13(1), 39-66. doi: https://doi.org/10.1080/15326349708807412
