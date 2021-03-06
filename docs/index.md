# BioCCP
## Intro
 
BioCCP.jl applies the Coupon Collector Problem to **combinatorial biotechnology**, in particular to aid (expected) **minimum sample size** determination for screening experiments. 

**Modular designs** are considered, created by randomly combining `r` modules from a set of `n` available modules (sampling with replacement). The module probabilities during the generation of the designs are specified by a probability/abundance vector `p`. Depending on how many complete sets of modules one wants to observe, parameter `m` can be increased from its default value of 1 to a higher value. 

For a specific combinatorial design set-up of interest, a report with results regarding (expected) minimum sample sizes can be easily retrieved by using the provided Interactive Pluto notebook. Additionally, a Pluto notebook with case studies is provided to illustrate the usage of BioCCP onto real biological examples.

>  **References:** 
>  
> Boneh, A., &amp; Hofri, M. (1997). The coupon-collector problem revisited—a survey of engineering problems and computational methods. *Stochastic Models, 13*(1), 39-66.
> 
> Doumas, A. V., &amp; Papanicolaou, V. G. (2016). The coupon collector’s problem revisited: generalizing the double Dixie cup problem of Newman and Shepp. *ESAIM: Probability and Statistics, 20*, 367-399.

## Functions

<html lang="en"><head><meta charset="UTF-8"/><meta name="viewport" content="width=device-width, initial-scale=1.0"/></head><body><div id="documenter"><nav class="docs-sidebar">  <article class="docstring"><header><a class="docstring-binding" id="BioCCP.expectation_minsamplesize" href="#BioCCP.expectation_minsamplesize"><code><strong>BioCCP.expectation_minsamplesize</strong></code></a> — <span class="docstring-category">Function</span></header><section><div><pre><code class="language-julia">expectation_minsamplesize(n::Integer; p=ones(n)/n, m::Integer=1, r=1, normalize=true)</code></pre><p>Calculates the expected minimum number of designs <code>E[T]</code>  to observe each module at least <code>m</code> times.</p><ul><li><code>n</code>: number of modules in the design space</li><li><code>p</code>: vector with the probabilities or abundances of the different modules</li><li><code>m</code>: number of times each module has to be observed in the sampled set of designs </li><li><code>r</code>: number of modules per design</li><li>normalize: if true, normalize <code>p</code></li></ul>
<pre><code class="language-julia-repl">julia&gt; n = 100
julia&gt; expectation_minsamplesize(n; p = ones(n)/n, m = 1, r = 1, normalize = true)
518</code></pre></div></section>



</article><article class="docstring">
 <header><a class="docstring-binding" id="BioCCP.std_minsamplesize" href="#BioCCP.std_minsamplesize"><code><strong>BioCCP.std_minsamplesize</strong></code></a> — <span class="docstring-category">Function</span></header><section><div><pre><code class="language-julia">std_minsamplesize(n::Integer; p=ones(n)/n, m::Integer=1, r=1, normalize=true)</code></pre><p>Calculates the standard deviation on the minimum number of designs to observe each module at least <code>m</code> times.</p><ul><li><code>n</code>: number of modules in the design space</li><li><code>p</code>: vector with the probabilities or abundances of the different modules</li><li><code>m</code>: number of times each module has to be observed in the sampled set of designs </li><li><code>r</code>: number of modules per design</li><li>normalize: if true, normalize <code>p</code></li></ul>
 
 <pre><code class="language-julia-repl">julia&gt; n = 100
julia&gt; std_minsamplesize(n; p = ones(n)/n, m = 1, r = 1, normalize = true)
126</code></pre></div></section>



</article><article class="docstring"><header><a class="docstring-binding" id="BioCCP.success_probability" href="#BioCCP.success_probability"><code><strong>BioCCP.success_probability</strong></code></a> — <span class="docstring-category">Function</span></header><section><div><pre><code class="language-julia">success_probability(n::Integer, t::Integer; p=ones(n)/n, m::Integer=1, r=1, normalize=true)</code></pre><p>Calculates the success probability <code>F(t) = P(T ≤ t)</code> or the probability that the minimum number of designs <code>T</code> to see each module at least <code>m</code> times is smaller than or equal to <code>t</code>.</p><ul><li><code>n</code>: number of modules in design space</li><li><code>t</code>: sample size or number of designs for which to calculate the success probability </li><li><code>p</code>: vector with the probabilities or abundances of the different modules</li><li><code>m</code>: number of times each module has to be observed in the sampled set of designs </li><li><code>r</code>: number of modules per design</li><li>normalize: if true, normalize <code>p</code></li></ul>
 <pre><code class="language-julia-repl">julia&gt; n = 100
julia&gt; t = 600
julia&gt; success_probability(n, t; p = ones(n)/n, m = 1, r = 1, normalize = true)
0.7802171997092149</code></pre></div></section>
 
 </article><article class="docstring"><header><a class="docstring-binding" id="BioCCP.expectation_fraction_collected" href="#BioCCP.expectation_fraction_collected"><code><strong>BioCCP.expectation_fraction_collected</strong></code></a> — <span class="docstring-category">Function</span></header><section><div><pre><code class="language-julia">expectation_fraction_collected(n::Integer, t::Integer; p=ones(n)/n, r=1, normalize=true)</code></pre><p>Calculates the fraction of the total number of modules that is expected to be observed after collecting <code>t</code> designs.</p><ul><li><code>n</code>: number of modules in design space</li><li><code>t</code>: sample size or number of designs </li><li><code>p</code>: vector with the probabilities or abundances of the different modules </li><li><code>r</code>: number of modules per design</li><li>normalize: if true, normalize <code>p</code></li></ul>
 
<pre><code class="language-julia-repl">julia&gt; n = 100
julia&gt; t = 200
julia&gt; expectation_fraction_collected(n, t; p = ones(n)/n, r = 1, normalize=true)
0.8660203251420364</code></pre></div></section>
 
 
 </article><article class="docstring"><header><a class="docstring-binding" id="BioCCP.prob_occurrence_module" href="#BioCCP.prob_occurrence_module"><code><strong>BioCCP.prob_occurrence_module</strong></code></a> — <span class="docstring-category">Function</span></header><section><div><pre><code class="language-julia">prob_occurrence_module(p<sub>i</sub>, t::Integer, r, k::Integer)</code>   </pre><p>Calculates the probability that specific module with module probability <code>p<sub>i</sub></code>  has occurred <code>k</code> times after collecting <code>t</code> designs.</p><p>Sampling processes of individual modules are assumed to be independent Poisson processes.</p><ul><li><code>p<sub>i</sub></code>: module probability</li><li><code>t</code>: sample size or number of designs </li><li><code>r</code>: number of modules per design </li><li><code>k</code>: the number of occurrences </li></ul>
 
<pre><code class="language-julia-repl">julia&gt; p<sub>i</sub> = 0.005
julia&gt; t = 500
julia&gt; r = 1
julia&gt; k = 2
julia&gt; prob_occurrence_module(p<sub>i</sub>, t, r, k)
0.25651562069968376</code></pre> 


 


