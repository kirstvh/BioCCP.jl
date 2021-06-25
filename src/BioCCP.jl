module BioCCP

#using 

# export all functions that are relevant for the user
export expectation_minsamplesize, std_minsamplesize, success_probability, 
        expectation_fraction_collected, prob_occurrence_module


"""
    exp_ccdf(n, T; p_vec=ones(n), m=1, r=1, normalize=true)

Calculates `1 - F(t)`, which is the complement of the success probability
`F(t)=P(T ≤ t)` (= probability that the expected minimum
  number of designs T is smaller than `t` in order to 
  see each module at least `m` times). This function
  serves as the integrand for calculating `E[T]`.
 
- `n`: number of modules in the design space
- `p_vec`: vector with the probabilities/abundances of the different modules in the design space during library generation
- `T`: number of designs
- `m`: number of times each module has to observed in the sampled set of designs
- `r`: number of modules per design
- normalize: if true, normalize `p_vec`

References:
- Doumas, A. V., & Papanicolaou, V. G. (2016). The coupon collector’s problem revisited:
        generalizing the double Dixie cup problem of Newman and Shepp. ESAIM: Probability and Statistics, 20, 367-399.
- Boneh, A., & Hofri, M. (1997). The coupon-collector problem revisited—a survey of
        engineering problems and computational methods. Stochastic Models, 13(1), 39-66.

## Examples
 

```julia-repl
julia> n = 100
julia> t = 500
julia> exp_ccdf(n, t; p_vec=ones(n), m=1, r=1, normalize=true)
0.4913906004535237
```
"""
function exp_ccdf(n, t; p_vec=ones(n), m=1, r=1, normalize=true)   
    @assert length(p_vec) == n
	    
    # Normalize probabilities
    if normalize
        p_vec=p_vec ./ sum(p_vec)    
    end   
    # Initialize probability P
    P_cdf = 1
    for i in 1:n
          Sm = 0
        for j in 1:m
            Sm += ((p_vec[i]*r*t)^(j-1))/factorial(j-1) #formulas see paper <Introduction
        end 
        P_cdf *= (1 - Sm*exp(-p_vec[i]*r*t))        
    end   
    P = 1 - P_cdf
    return P
end   

"""
    approximate_moment(n, fun; p_vec=ones(n), q=1, m=1, r=1,
steps=10000, normalize=true)

Calculates the q-th rising moment of `T[N]` (number of designs that are needed to collect
all modules `m` times). Integral is approximated by the Riemann sum.

Reference: 
- Doumas, A. V., & Papanicolaou, V. G. (2016). The coupon collector’s problem revisited:
        generalizing the double Dixie cup problem of Newman and Shepp. ESAIM: Probability
        and Statistics, 20, 367-399.

## Examples

```julia-repl
julia> n = 100
julia> fun = exp_ccdf
julia> approximate_moment(n, fun; p_vec=ones(n), q=1, m=1, r=1,
steps=10000, normalize=true)
518.8175339489885
```
"""
function approximate_moment(n, fun; p_vec=ones(n)/n, q=1, m=1, r=1,
	        steps=10000, normalize=true)
    @assert length(p_vec) == n
    a = 0; b = 0
    while fun(n, b; p_vec=p_vec, m=m, r=r, normalize=normalize) > 0.00001
        b += 5
    end
    δ = (b-a)/steps; t = a:δ:b
    #integration exp_ccdf, see paper References [1]
    qth_moment = q .* sum(δ .* fun.(n, t; p_vec=p_vec, m=m, r=r, normalize=normalize) .* t.^[q-1]) 
    return qth_moment           
end

"""
    expectation_minsamplesize(n; p_vec=ones(n), m=1, r=1, normalize=true)

Calculates the expected minimum number of designs  `E[T]` to observe each module at least `m` times.

- `n`: number of modules in the design space
- `p_vec`: vector with the probabilities or abundances of the different modules
- `m`: number of times each module has to be observed in the sampled set of designs 
- `r`: number of modules per design
- normalize: if true, normalize `p_vec`

References:
- Doumas, A. V., & Papanicolaou, V. G. (2016). The coupon collector’s problem revisited:
    generalizing the double Dixie cup problem of Newman and Shepp. ESAIM: Probability
    and Statistics, 20, 367-399.
- Boneh, A., & Hofri, M. (1997). The coupon-collector problem revisited—a survey of
    engineering problems and computational methods. Stochastic Models, 13(1), 39-66.


## Examples

```julia-repl
julia> n = 100
julia> expectation_minsamplesize(n; p_vec=ones(n), m=1, r=1, normalize=true)
519.0
```
"""
function expectation_minsamplesize(n; p_vec=ones(n), m=1, r=1, normalize=true)
    @assert length(p_vec) == n
    E = approximate_moment(n, exp_ccdf; p_vec=p_vec, q=1, m=m, r=r, normalize=normalize)
    return ceil(E)
end

"""
    std_minsamplesize(n; p_vec=ones(n), m=1, r=1, normalize=true)

Calculates the standard deviation on the minimum number of designs to observe each module at least `m` times.
    
- `n`: number of modules in the design space
- `p_vec`: vector with the probabilities or abundances of the different modules
- `m`: number of complete sets of modules that need to be collected 
- `r`: number of modules per design
- normalize: if true, normalize `p_vec`

## Examples

```julia-repl
julia> n = 100
julia> std_minsamplesize(n; p_vec=ones(n), m=1, r=1, normalize=true)
126.0
```
"""
function std_minsamplesize(n; p_vec=ones(n), m=1, r=1, normalize=true)
    @assert length(p_vec) == n
    M1 = approximate_moment(n, exp_ccdf; p_vec=p_vec, q=1, m=m, r=r,  normalize=normalize)
    M2 = approximate_moment(n, exp_ccdf; p_vec=p_vec, q=2, m=m, r=r, normalize=normalize)
    var = M2 - M1 - M1^2
    return ceil(sqrt(var))
end

"""
    success_probability(n, t; p_vec=ones(n), m=1, r=1, normalize=true)

Calculates the success probability `F(t) = P(T ≤ t)` or the probability that 
the minimum number of designs `T` to see each module at least `m` times
is smaller than `t`.

- `n`: number of modules in design space
- `t`: sample size/number of designs for which to calculate the success probability 
- `p_vec`: vector with the probabilities or abundances of the different modules
- `m`: number of complete sets of modules that need to be collected 
- `r`: number of modules per design
- normalize: if true, normalize `p_vec`

References:
- Boneh, A., & Hofri, M. (1997). The coupon-collector problem revisited—a survey of engineering problems and computational methods. Stochastic Models, 13(1), 39-66.

## Examples

```julia-repl
julia> n = 100
julia> t = 600
julia> success_probability(n, t; p_vec=ones(n), m=1, r=1, normalize=true)
0.7802171997092149
```
"""
function success_probability(n, t; p_vec=ones(n), m=1, r=1, normalize=true)   
    P_success = 1 - exp_ccdf(n, t; p_vec=p_vec, m=m, r=r, normalize=normalize) 
    return P_success
end

"""
    expectation_fraction_collected(n, t; p_vec=ones(n), r=1, normalize=true) 

Calculates the fraction of all modules that is expected to be observed
after collecting `t` designs.

- `n`: number of modules in design space
- `t`: sample size/number of designs 
- `p_vec`: vector with the probabilities or abundances of the different modules 
- `r`: number of modules per design
- normalize: if true, normalize `p_vec`

References:
- Boneh, A., & Hofri, M. (1997). The coupon-collector problem revisited—a survey of engineering problems and computational methods. Stochastic Models, 13(1), 39-66.

## Examples

```julia-repl
julia> n = 100
julia> t = 200
julia> expectation_fraction_collected(n, t; p_vec=ones(n), r=1, normalize=true)
0.8660203251420364
```
"""
function expectation_fraction_collected(n, t; p_vec=ones(n), r=1, normalize=true)
    if normalize
        p_vec = p_vec./sum(p_vec)
    end
    frac = sum( (1-(1-p_vec[i])^(t*r)) for i in 1:n )/n
    return frac
end

"""
    prob_occurrence_module(p, t, j)

Calculates probability that specific module with module probability `p` 
has occurred `k` times after collecting `t` designs.

Sampling processes of individual modules are assumed to be independent Poisson processes.

- `p`: module probaility
- `t`: sample size/number of designs 
- `k`: number of occurrence 

References:
- Boneh, A., & Hofri, M. (1997). The coupon-collector problem revisited—a survey of engineering problems and computational methods. Stochastic Models, 13(1), 39-66.

## Examples

```julia-repl
julia> p = 0.005
julia> t = 500
julia> k = 2
julia> prob_occurrence_module(p, t, k)
0.25651562069968376
```
"""
function prob_occurrence_module(p, t, k)
	return (exp(-1*(p*t))*(p*t)^k)/factorial(k) 
end

end