### A Pluto.jl notebook ###

# v0.12.15

using Markdown
using InteractiveUtils

# ╔═╡ 5ad5c202-20f8-11eb-23f1-4f38b687c285
using STMOZOO.Example

# ╔═╡ d2d007b8-20f8-11eb-0ddd-1181d4565a85
using Plots

# ╔═╡ 45189a82-20fa-11eb-0423-05ce1b84639d
using Zygote

# ╔═╡ 171fee18-20f6-11eb-37e5-2d04caea8c35
md"""
# Example: solving quadratic systems

By Michiel Stock

## Introduction

In this notebook, I have chosen to write code so minimize *quadratic functions* of the form:

$$f(\mathbf{x}) =  \frac{1}{2} \mathbf{x}^\intercal P\mathbf{x} + \mathbf{q} \cdot \mathbf{x} + r\,,$$

where $P \succ 0$.

Hence, solving

$$\min_\mathbf{x}\, f(\mathbf{x})\,.$$
"""

# ╔═╡ fb4aeb8c-20f7-11eb-0444-259de7b76883
md"""

## Derivation

The gradient of the quadratic function is

$$\nabla f(\mathbf{x})=P\mathbf{x} +\mathbf{q}\,.$$

Setting this to zero gives

$$\mathbf{x}^\star=-P^{-1}\mathbf{q}\,.$$

We also know that $\nabla^2f(\mathbf{x}) = P \succ 0$, so $\mathbf{x}^\star$ is a minimizer.
"""

# ╔═╡ 52e87238-20f8-11eb-2ea0-27ee1208d3c3
md"""
## Illustration

We will use a simple 2D example to illustrate this code.
"""

# ╔═╡ e9a9d69e-20f8-11eb-0d9d-330ee5e9cf25
P = [7 -1; -1 2]

# ╔═╡ fdd4e550-20f8-11eb-227b-25f36708484d
q = [-14.0, -2.0]

# ╔═╡ 025fd6e8-20f9-11eb-3e7d-3519f3c4b58f
r = 1

# ╔═╡ 096eff98-20f9-11eb-1e61-99d5714895ba
md"We have a function to generate the quadratic function."

# ╔═╡ 165509ca-20f9-11eb-107c-550cbba0f0e9
f_quadr = quadratic_function(P, q, r)

# ╔═╡ 1fffc82a-20f9-11eb-198c-c160d7dac87d
f_quadr([2, 1])

# ╔═╡ 26ab6ce2-20f9-11eb-1836-1756b290e5e3
md"No more need to remember the formulla for the minimizer! Just use `solve_quadratic_system`!"

# ╔═╡ 49832a8e-20f9-11eb-0841-19a40a12db18
x_star = solve_quadratic_system(P, q, r)

# ╔═╡ 55e0e274-20f9-11eb-36c0-753f228f7e9b
begin
	contourf(-20:0.1:20, -20:0.1:20, (x, y) -> f_quadr([x,y]), color=:speed)
	scatter!([x_star[1]], [x_star[2]], label="minimizer")
end

# ╔═╡ b1551758-20f9-11eb-3e8f-ff9a7127d7f8
md"""
## Approximating non-quadratic functions

We can approximate non-quadratic functions by a quadratic function: The second order Taylor approximation $\hat{f}$ of a function $f$ at $\mathbf{x}$ is
$$f(\mathbf{x}+\mathbf{v})\approx\hat{f}(\mathbf{x}+\mathbf{v}) = f(\mathbf{x}) + \nabla f(\mathbf{x})^\top \mathbf{v} + \frac{1}{2} \mathbf{v}^\top \nabla^2 f(\mathbf{x}) \mathbf{v}\,.$$

Let us use this idea for the Rosenbrock function.
"""

# ╔═╡ 41d8f1dc-20fa-11eb-3586-a989427c1fd6
f_nonquadr((x1, x2); a=1, b=5) = (a-x1)^2 + b * (x2 - x1^2)^2

# ╔═╡ 4ed4215e-20fa-11eb-11ee-f7741591163c
x = [0.0, 0.0]

# ╔═╡ 56af99ee-20fa-11eb-0240-69c675efb78c
fx = f_nonquadr(x)

# ╔═╡ 6c5473b4-20fa-11eb-327b-51ac560530eb
∇fx = f_nonquadr'(x)

# ╔═╡ 7518c2c0-20fa-11eb-32c0-a9db2a91cbc5
∇²fx = Zygote.hessian(f_nonquadr, x)

# ╔═╡ 34027942-20fb-11eb-261e-3b991ce4c9f8
v = solve_quadratic_system(∇²fx, ∇fx, fx)

# ╔═╡ 3bbeb85c-20fc-11eb-04d0-fb12d8ace50a
f̂(x′) = quadratic_function(∇²fx, ∇fx, fx)(x′ .- x)

# ╔═╡ 8623ac1a-20fa-11eb-2d45-49cce0fdac86
begin
	plot_nonquadr = contourf(-2:0.01:2, -2:0.01:2, (x, y) -> f_nonquadr([x,y]), color=:speed, title="non-quadratic function")
	scatter!(plot_nonquadr, [x[1]], [x[2]], label="x")
	scatter!(plot_nonquadr, [x[1]+v[1]], [x[2]+v[2]], label="x + v")
	
	plot_approx = contourf(-2:0.01:2, -2:0.01:2, (x, y) -> f̂([x,y]), color=:speed,
		title="quadratic approximation")
	scatter!(plot_approx, [x[1]], [x[2]], label="x")
	scatter!(plot_approx, [x[1]+v[1]], [x[2]+v[2]], label="x + v")
	
	plot(plot_nonquadr, plot_approx, layout=(2,1), size=(600, 800))
	
end


# ╔═╡ Cell order:
# ╟─171fee18-20f6-11eb-37e5-2d04caea8c35
# ╠═5ad5c202-20f8-11eb-23f1-4f38b687c285
# ╠═d2d007b8-20f8-11eb-0ddd-1181d4565a85
# ╟─fb4aeb8c-20f7-11eb-0444-259de7b76883
# ╟─52e87238-20f8-11eb-2ea0-27ee1208d3c3
# ╠═e9a9d69e-20f8-11eb-0d9d-330ee5e9cf25
# ╠═fdd4e550-20f8-11eb-227b-25f36708484d
# ╠═025fd6e8-20f9-11eb-3e7d-3519f3c4b58f
# ╟─096eff98-20f9-11eb-1e61-99d5714895ba
# ╠═165509ca-20f9-11eb-107c-550cbba0f0e9
# ╠═1fffc82a-20f9-11eb-198c-c160d7dac87d
# ╟─26ab6ce2-20f9-11eb-1836-1756b290e5e3
# ╠═49832a8e-20f9-11eb-0841-19a40a12db18
# ╠═55e0e274-20f9-11eb-36c0-753f228f7e9b
# ╟─b1551758-20f9-11eb-3e8f-ff9a7127d7f8
# ╠═41d8f1dc-20fa-11eb-3586-a989427c1fd6
# ╠═45189a82-20fa-11eb-0423-05ce1b84639d
# ╠═4ed4215e-20fa-11eb-11ee-f7741591163c
# ╠═56af99ee-20fa-11eb-0240-69c675efb78c
# ╠═6c5473b4-20fa-11eb-327b-51ac560530eb
# ╠═7518c2c0-20fa-11eb-32c0-a9db2a91cbc5
# ╠═34027942-20fb-11eb-261e-3b991ce4c9f8
# ╠═3bbeb85c-20fc-11eb-04d0-fb12d8ace50a
# ╟─8623ac1a-20fa-11eb-2d45-49cce0fdac86
