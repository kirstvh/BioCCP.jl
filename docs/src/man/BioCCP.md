# BioCCP
## Introduction 
 
BioCCP.jl applies the Coupon Collector's Problem to combinatorial biotechnology, in particular to help determining minimum sample sizes of screening experiments. Modular designs are considered, created by randomly combining `r` modules from a set of `n`available modules. The module proabbilities during the generation of the designs are specified by a probability/abundance vector `p_vec`. Depending on how many complete sets of modules one wants to observe, parameter `m` can be increased from its default value of 1 to the desired value. 

For a specific combinatorial design set-up of interest, a report with results regarding minimum sample sizes can be easily retrieved by using the provided Pluto notebook.

## Functions
```@docs
expectation_minsamplesize
std_minsamplesize
success_probability
expectation_fraction_collected
prob_occurence_module
```


## References
- Doumas, A. V., & Papanicolaou, V. G. (2016). The coupon collector’s problem revisited: generalizing the double Dixie cup problem of Newman and Shepp. ESAIM: Probability and Statistics, 20, 367-399.
- Boneh, A., & Hofri, M. (1997). The coupon-collector problem revisited—a survey of engineering problems and computational methods. Stochastic Models, 13(1), 39-66.




