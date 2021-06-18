# BioCCP
*Collecting Coupons in combinatorial biotechnology*

## Intro
During the **combinatorial engineering of biosystems** such as proteins, genetic circuits and genomes, diverse libraries are generated by assembling and recombining modules. The variants with the optimal phenotype are selected with screening techniques. However, when the number of available modules to compose a biological designs increases, a combinatorial explosion of design possibilities arises, allowing only for a part of the libary to be analyzed. In this case, it is important for a researcher to get insight in which (minimum) sample size sufficiently covers the design space, i.e. what is the expected **minimum number of designs so that all modules are observed at least once**.


## Functions
BioCCP contains functions for calculating minimum sample sizes and related statistics:

Function name    | Short description
---------------- | -----------------
`expectation_minsamplesize`        | Calculates the expected minimum number of designs to observe all modules  
`std_minsamplesize`      | Calculates standard deviation on minimum number of designs 
`success_probability`         | Calculates the probability that the minimum number of designs T is smaller than or equal to a given sample size t  
`expectation_fraction_collected` | Returns for a given sample size the fraction of the available modules that is expected to be observed
`prob_occurrence_module` | Calculates for a module with specified module probability *p*, the probability that the module has occurred *j* times when a given number of designs has been collected
 

For more info about the implementation of the functions, please consult the [docs](https://kirstvh.github.io/BioCCP/).

## Pluto notebook

The Pluto notebook provides an interactive illustration of all functions in BioCCP and assembles a report for your specific design set-up. 

Inputs for generating the report    |  
---------------- | 

Symbol    | Short description
---------------- | -----------------
 *n*       |  Total number of modules in the design space
*r*     |  The number of modules per design 
 *m*        | The number of times each module has to be observed (default = 1)
 *p<sub>vec</sub>*   (\*) |  Probability distribution of the modules 

>  (\*) 
>  *When exact probabilities are known*, define your custom module probability/abundance vector or load them in the notebook from an external file.
>  *When probabilities and/or their distribution are unknown*, you can either:
 > >  1) Assume the probabilities of all modules to be equal (uniform distribution), or
 > >  2) Assume the module probabilities to follow *Zipf's law*, specifying p<sub>max</sub> and p<sub>min</sub>, or
 > >  3) Assume the histogram of the module probabilities to behave like a *bell curve*, specifying the mean and variance       

Using the inputs, a report for sample size determination is created using the [functions](https://kirstvh.github.io/BioCCP/) described above. The report contains the following sections:

Report section    |   Short description       
---------------- |  -----------------           
Module probabilities       |     This section shows a plot with the probability of each module in the design space during library generation.   
Expected minimum sample size      |     This part displays the expected minimum number of designs *E[T]* and the standard deviation *std[T]*.         
Success probability      |    In this section, the report calculates the probability *F(t)* that the minimum number of designs *T* is smaller than or equal to a given sample size *t*. Moreover, a curve describing the success probability *F(t)* in function of an increasing sample size *t* is available, to determine a minimum sample size according to a probability cut-off.  
Expected observed fraction of the total number of modules        |    Here, the fraction of the total number of modules in the design space that is expected to be observed is computed for a given sample size *t*. A saturation curve, displaying the expected fraction of modules observed in function of increasing sample size, is provided.
Number of occurrence of a specific module      |      In this last part, you can specify the probability *p<sub>j</sub>* of the module of interest together with a particular sample size, to calculate a curve showing the probability for a module to occur *k* times (in function of *k*).   


## Getting started

#### Launch Pluto notebook from Browser 

Launch the Pluto notebook directly from your browser using Binder (no installation of Julia/packages required): [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kirstvh/PlutoNotebooks/main?urlpath=pluto/open?path=/home/jovyan/notebooks/BioCCP_Interactive_Notebook.jl)

#### Execute functions in Julia

&emsp; **(1)** &emsp; [Install Julia](https://julialang.org/downloads/) 

&emsp; **(2)** &emsp; Install BioCCP in the Julia REPL (hit ] to enter package mode):

    ] add https://github.com/kirstvh/BioCCP

&emsp; **(3)** &emsp; For using the [Pluto notebook](BioCCP/notebooks/BioCCP_Interactive_Notebook.jl):

&emsp;&emsp; In the Julia REPL, hit the following command to install the additional packages [Pluto](https://github.com/fonsp/Pluto.jl), [PlutoUI](https://github.com/fonsp/PlutoUI.jl) and [Plots](https://github.com/JuliaPlots/Plots.jl) 
  
    ] add Pluto, PlutoUI, Plots

&emsp;&emsp; Then start Pluto in the Julia REPL:

    using Pluto; Pluto.run()
    
&emsp;&emsp; Open the [notebook file](/notebooks/BioCCP_Interactive_Notebook.jl).

## References
Implementation of formula's was based on:

> Doumas, A. V., & Papanicolaou, V. G. (2016). *The coupon collector’s problem revisited: generalizing the double Dixie cup problem of Newman and Shepp.* ESAIM: Probability and Statistics, 20, 367-399. doi: 	https://doi.org/10.1051/ps/2016016

> Boneh, A., & Hofri, M. (1997). *The coupon-collector problem revisited—a survey of engineering problems and computational methods.* Stochastic Models, 13(1), 39-66. doi: https://doi.org/10.1080/15326349708807412
