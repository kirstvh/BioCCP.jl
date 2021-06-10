# Michiel Stock
# Example of a source code file implementing a module.


# all your code is part of the module you are implementing
module Example

# you have to import everything you need for your module to work
# if you use a new package, don't forget to add it in the package manager
using LinearAlgebra

# export all functions that are relevant for the user
export solve_quadratic_system, quadratic_function

"""
    solve_quadratic_system(P, q, r=0; testPD=false)

Find the minimizer of a canonical quadratic system:

``f(x) = 0.5xᵀ  P  x + x ⋅ q + r``

with `P` a positive scalar or a positive-definite n x n matrix and `q` a scalar or
a n x n vector. The intercept `r` can optionally be given but does not influence
the minimizer. One can use the keyword argument `testPD` to test if `P` is positive definite,
though this comes at a computational cost.

## Examples

Scalar case:

```julia-repl
julia> solve_quadratic_system(8, -4, 3)
0.5
```

Vector case:

```julia-repl
julia> P = [3 1; 1 2];

julia> q = [0.5, -2];

julia> solve_quadratic_system(P, q)
2-element Array{Float64,1}:
 -0.6               
  1.2999999999999998
```
"""
function solve_quadratic_system(P, q, r=0; testPD=false)
    testPD && @assert isposdef(P) "Provided P is not positive definite"
    return - P \ q
end

"""
    quadratic_function(P, q, r)

Give the paramters of a canonical quatratic function and returns a function.
"""
quadratic_function(P, q, r) = x -> 0.5x' * P * x + q ⋅ x + r 

end