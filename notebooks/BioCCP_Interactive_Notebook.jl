### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 2d3ad982-ef1f-45ae-b247-9679c0faa853
using BioCCP, Plots, PlutoUI

# ╔═╡ 4d246460-af05-11eb-382b-590e60ba61f5
md"## Collecting Coupons in combinatorial biotechnology

This notebook provides functions and visualizations to determine expected minimum sample sizes for biotechnological experiments based on the mathematical framework of the Coupon Collector Problem (implemented formulas based on [^1], [^2]).

"

# ╔═╡ a2fd6000-1450-4dfe-9426-5303ae64bfb3
md"""Please install the packages `BioCCP`, `Plots` and `PlutoUI` in the Julia Package Manager for this notebook to work."""

# ╔═╡ a8c81622-194a-443a-891b-bfbabffccff1
begin
md""" 
 
👇 **COMPLETE THE FIELDS BELOW** 👇

№ modules in design space:                       $(@bind n_string TextField(default = "100")) \
	
№ modules per design:                            $(@bind r NumberField(1:20))\
№ complete sets of modules to collect:               $(@bind m NumberField(1:20))\
	
Abundances of modules during library generation:       $(@bind ps Select(["Equal", "Unequal"], default = "Equal"))"""
	
end

# ╔═╡ 45507d48-d75d-41c9-a018-299e209f900e
begin
	n = parse(Int64, n_string);
	if ps == "Equal"
		distribution = "Equal"
	end
		if ps == "Unequal"	
	md""" 	                         ↳     Specify distribution:                         
	$(@bind distribution Select(["Bell curve", "Zipf's law", "Custom vector"], default = " "))"""
		end	
end

# bell curve ipv normale distrbution, neem quantielen, niet samplen, vaste uitkomst

# ╔═╡ b17f3b8a-61ee-4563-97cd-19ff049a8e1e
begin
	if distribution == "Bell curve"					
			md"""                                          pₘₐₓ/pₘᵢₙ:  $(@bind pmaxpmin_str TextField(default = "4")) 
                                            """
			end
end

# ╔═╡ e3b4c2d8-b78c-467e-a863-5eecb8ec58dc
begin
	if distribution == "Zipf's law"
		md"""                                       pₘₐₓ/pₘᵢₙ:   $(@bind pmaxpmin_string TextField(default = "4")) 
                                            """
			end

end

# ╔═╡ 2639e3fb-ccbb-44de-bd15-1c5dbf6c1539
begin
	if distribution == "Custom vector"
				md"""             	       ↳  Enter/load your custom abundances by changing the cell below 👇"""			
		end
end

# ╔═╡ 44d4dfee-3073-49aa-867c-3abea10e6e37
begin
	# load your custom probability vector!
	# see for example XLSX package
	if distribution == "Custom vector"
		abundances = rand(200:1:400, n)
	end
end

# ╔═╡ f6ebf9fb-0a29-4cb4-a544-6c6e32bedcc4
md"""	
 
🎯 **REPORT**  🎯

**💻 Module probabilities**                                                                                                                       $(@bind show_modprobs Select(["🔻 SHOW ", "🔺 HIDE "], default="🔺 HIDE ") )  \
*How the abundances of the modules are distributed during combinatorial library generation.*
"""

# ╔═╡ d4a9da7a-f455-426b-aecd-227c25e1d4e8
begin

	if ps == "Equal"
	 	
		p_vec = ones(n)./sum(ones(n));
		
	elseif ps == "Unequal"
		if distribution == "Bell curve"
			ratio = parse(Float64, pmaxpmin_str)
			ab1 = 1
			ab2 = ratio*ab1
			μ = (ab1+ab2)/2
			σ = (ab2-ab1)/6
			
			#create fixed distribution of abundances according to percentiles of bell curve
			n_perc_1 = Int(floor(n*0.34)); 
			n_perc_2 = Int(floor(n*0.135));
			n_perc_3 = Int(floor(n*0.0215));
			#n_perc_4 = Int(floor(n*0.0013));
			n_perc_rest = n - 2*n_perc_1 - 2*n_perc_2 - 2*n_perc_3 ;
			p_vec_unnorm = vcat(fill(μ,2*n_perc_1+n_perc_rest), fill(μ+1.5*σ, n_perc_2), fill(μ-1.5*σ, n_perc_2), fill(μ+3*σ, n_perc_3), fill(μ-3*σ, n_perc_3) )
		
			# normalize sum to 1
			p_vec = sort(p_vec_unnorm ./ sum(p_vec_unnorm))
		end
		
		if distribution == "Custom vector"
			p_vec_unnorm = abundances
			p_vec = abundances ./ sum(abundances)
		end
		
		if distribution == "Zipf's law"
			ratio = parse(Float64, pmaxpmin_string)
			α = exp(log(ratio)/(n-1))
			p_vec = collect(α.^-(1:n))
			p_vec = p_vec ./ sum(p_vec)
		end
	end
	
	if show_modprobs == "🔻 SHOW "   
	
	scatter(p_vec, title = "Probability mass function", ylabel = "module probability pⱼ", xlabel = "module j", label="", size = (700, 400))
	ylims!((0,maximum(p_vec) + maximum(p_vec)-minimum(p_vec) ))

	end	
end

# ╔═╡ 87c3f5cd-79bf-4ad8-b7f8-3e98ec548a9f
begin
	if show_modprobs == "🔻 SHOW "  && distribution == "Bell curve"
		histogram(p_vec, normalize=:probability,  bar_edges=false,  size = (650, 340), orientation=:v, bins=[(μ -  3.2*σ)/sum(p_vec_unnorm), (μ - 2*σ)/sum(p_vec_unnorm), (μ-σ)/sum(p_vec_unnorm), (μ + σ)/sum(p_vec_unnorm), (μ + 2*σ)/sum(p_vec_unnorm), (μ +  3.2*σ)/sum(p_vec_unnorm)])
		# if distribution == "Normally distributed"
		# 	plot!(x->pdf(Normal(μ, σ), x), xlim=xlims())
		# 	xlabel!("Abundance"); ylabel!("probability"); title!("Distribution of module abundances")
		# end
		xlabel!("Probability"); ylabel!("Relative frequency"); title!("Distribution of module probabilities")
	end	
end

# ╔═╡ d877bd4c-497d-46d1-9c58-b6fe26933bfc
begin
	if show_modprobs == "🔻 SHOW "  && distribution == "Bell curve"
md"""For $n_string modules of which the probabilities form a bell curve with ratio pₘₐₓ/pₘᵢₙ = $pmaxpmin_str , we follow the percentiles of a normal distribution to generate the probability vector.

We consider μ to be the mean module probability and σ to be the standard deviation of the module probabilities.
		
According to the percentiles
- 68% of the module probabilities lies in the interval [μ - σ, μ + σ], 
- 95% of falls into the range [μ - 2σ, μ + 2σ] and 
- 99.7% lies in [μ - 3σ, μ +3σ]. 
		
We use the ratio pₘₐₓ/pₘᵢₙ to fix the width of the interval [μ - 3σ, μ +3σ]. (We assume that pₘₐₓ = μ +3σ and pₘᵢₙ = μ - 3σ and calculate μ and σ from this assumption). In addition, we make sure the sum of the probability vector sums up to 1.
		
As a result, we get:
-  $(n_perc_1+n_perc_rest) modules with a probability of $(µ/sum(p_vec_unnorm))
-  $(n_perc_2)  modules with a probability of $((μ+1.5*σ)/sum(p_vec_unnorm))
-  $(n_perc_2)  modules with a probability of $((μ-1.5*σ)/sum(p_vec_unnorm))
-  $(n_perc_3)  modules with a probability of $((μ+2.5*σ)/sum(p_vec_unnorm))
-  $(n_perc_3)  modules with a probability of $((μ-2.5*σ)/sum(p_vec_unnorm))"""
	end	
end

# ╔═╡ f098570d-799b-47e2-b692-476a4d95825b
if show_modprobs == "🔻 SHOW " 
md"Each biological design in the design space is built by choosing $r module(s) (with replacement) out of a set of $n_string modules according to the module probabilities visualized above."
end

# ╔═╡ caf67b2f-cc2f-4d0d-b619-6e1969fabc1a
md""" **💻 Expected minimum sample size**                                                                                                             $(@bind show_E Select(["🔻 SHOW ", "🔺 HIDE "], default="🔺 SHOW ")) 
\
*The expected minimum number of designs to observe each module at least $m times in the sampled set of designs.* """  

# ╔═╡ 6f14a72c-51d3-4759-bb8b-10db1dc260f0
begin
	if show_E == "🔻 SHOW "   
		E = Int(expectation_minsamplesize(n; p_vec = p_vec, m=m, r = r))
		sd = Int(std_minsamplesize(n; p_vec = p_vec, m=m, r = r))
		
			md""" 
     `Expected minimum sample size`     = **$E designs**\
		
     `Standard deviation             `                           = **$sd designs**  	"""
	end
	# begin
		
	# 	#E_vec = []
	# 	#sd_vec = []
			 
	# 	#if ps == "Unequal" && probs_unequal_norm
	# 			#iter = 10
	# 			#for i in 1:iter
	# 				#p_vec_i = rand(Normal(μ,  σ), n)
	# 				#p_vec_i = p_vec_i ./ sum(p_vec_i)
	# 				#E_i = expectation_minsamplesize(n; p_vec = p_vec_i, m=m, q = q)
	# 				#sd_i = std_minsamplesize(n; p_vec = p_vec_i, m=m, q=q)
	# 				#push!(E_vec, E_i)
	# 				#push!(sd_vec, sd_i)
	# 			#end
	# 			#E = Int(ceil(mean(E_vec)))
	# 			#sd = Int(ceil(mean(sd_vec)))
	# 			#E_CI_lhs = Int(ceil( E - quantile(Normal(), 1-0.05/2)*sd/sqrt(iter)))
	# 			#E_CI_rhs = Int(ceil(E + quantile(Normal(), 1-0.05/2)*sd/sqrt(iter)))
	# 			#sd_CI_lhs = Int(ceil( sd - quantile(Normal(), 1-0.05/2)*sd/sqrt(iter)))
	# 			#sd_CI_rhs = Int(ceil(sd + quantile(Normal(), 1-0.05/2)*sd/sqrt(iter)))
				
			
	# 		#md""" 
	# 		#``` 
	# 		#Expected minimum sample size E[Tp]
	# 		#```		 	
	# 		#= **$E designs**    -------- 95% CI :  [$E_CI_lhs, $E_CI_rhs]
		
	# 		#``` 
	# 		#Standard deviation sd[Tp]  
	# 		#```	
		
	# 		#= **$sd designs**   -------- 95% CI :  [$sd_CI_lhs, $sd_CI_rhs]
		
	# 		#---------------
	# 		#"""
		
	# 		#else
	
	# #	E = expectation_minsamplesize(n; p_vec = p_vec, m=m, q = q)
	# #	sd = std_minsamplesize(n; p_vec = p_vec, m=m, q=q)
			
			
	# #		md""" 
	# #		``` 
	# #		Expected minimum sample size E[Tp]
	# #		```		 	
	# #		= **$E designs**    
		
	# #		``` 
	# #		Standard deviation sd[Tp]  
	# #		```	
		
	# #		= **$sd designs**   
		
	# #		---------------
	# #		"""
				 
	# 		#end
			
			
		 
	# #end
end

# ╔═╡ 3c07dd11-5be7-4ec7-992e-38dd07260d71


# ╔═╡ 22fe8006-0e81-4e0a-a460-28610a55cd97
md""" **💻 Success probability**                                                                                                                  $(@bind show_success Select(["🔻 SHOW ", "🔺 HIDE "], default="🔺 HIDE ") )\

*The probability that the minimum number of designs T is smaller than or equal to a given sample size t.* """

# ╔═╡ db4371e4-7f86-4db3-b076-12f6cd220b89
begin
	if show_success == "🔻 SHOW " 
		md"""    👉 Enter your sample size of interest: $(@bind sample_size_1_string TextField(default="500"))""" 
		
	end
	#genereer tabel + download knop
end

# ╔═╡ 317995ed-bdf4-4f78-bd66-a39ffd1dc452
begin
	if show_success == "🔻 SHOW " 
	sample_size_1 = parse(Int64, sample_size_1_string);
	
	p_success = success_probability(n, sample_size_1; p_vec = p_vec, m = m, r = r)
	
	md""" 
              ↳ `Success probability F(t)`  = **$p_success**\
	"""
	end
end

# ╔═╡ 9c36cfea-91dd-4bf7-b4f1-b34bce20cb9c


# ╔═╡ ca5a4cef-df67-4a5e-8a86-75a9fe8c6f37
if show_success == "🔻 SHOW " 
	md"*A curve describing the success probability in function of sample size.*"
end

# ╔═╡ 24f7aae7-d37a-4db5-ace0-c910b178da88
begin
if show_success == "🔻 SHOW " 
	
sample_size_initial = 5
	while (1 - success_probability(n, sample_size_initial; p_vec = p_vec, r = r, m = m)) > 0.0005
		global sample_size_initial += n/10
	end
		
	sample_sizes = 0: n/10 :sample_size_initial
	successes = success_probability.(n, sample_sizes; p_vec = p_vec, r = r, m = m)
plot(sample_sizes, successes, title = "Success probability in function of sample size", xlabel = "sample size s", ylabel= "P(s ≤ Sₘᵢₙ)", label = "", legend=:bottomright, size=(600,400), seriestype=:scatter )
		end
	 
end

# ╔═╡ 37f951ee-885c-4bbe-a05f-7c5e48ff4b6b
begin
	#following one-sided version of Chebyshev's inequality.
	 
	function chebyshev_onesided_larger(X, μ, σ)
		X_μ = X - μ
		return σ^2 / (σ^2 + X_μ^2)
	end
	function chebyshev_onesided_smaller(X, μ, σ)
		X_μ = μ - X
		return σ^2 / (σ^2 + X_μ^2)
	end
if show_success == "🔻 SHOW "
if sample_size_1 < E
	compare = "smaller"
		if sample_size_1 <= n/r
			print_sentence = "P(minimum sample size ≤ $sample_size_1) = 0."        
		else
	prob_chebyshev = chebyshev_onesided_smaller(sample_size_1, E, sd)
	print_sentence = "P(minimum sample size ≤ $sample_size_1) ≤ $prob_chebyshev. "
		end
		
elseif sample_size_1 > E
	compare = "greater"
	prob_chebyshev = chebyshev_onesided_larger(sample_size_1, E, sd)
	print_sentence = "P(minimum sample size ≥ $sample_size_1) ≤ $prob_chebyshev. "	
		
	elseif sample_size_1==E
		print_sentence = "P(minimum sample size ≤ $sample_size_1 OR minimum sample size ≥ $sample_size_1) ≤ 1."
		
end

	md"""*Upper bound on probability that minimum sample size is smaller than given sample size t, according to Chebychev's inequality*:
	                                                                                                                                            
		$print_sentence"""
	end
end

# ╔═╡ ca33610c-1be9-4c01-b0b7-ce4b2f7896df


# ╔═╡ dc696281-7a5b-4568-a4c2-8dde90af43f0
md""" **💻 Expected observed fraction of the total number of modules**                 $(@bind show_satur Select(["🔻 SHOW ", "🔺 HIDE "], default="🔺 HIDE "))\
*The expected fraction of the total number of available modules observed after collecting a given number of designs.*"""

# ╔═╡ eb92ff7c-0140-468c-8b32-f15d1cf15913
if show_satur == "🔻 SHOW " 
		md"""
   👉 Enter your sample size of interest: $(@bind sample_size_2_string TextField(default="50")) """ 
end

# ╔═╡ f0eaf96b-0bc0-4194-9a36-886cb1d66e00
begin
	if show_satur == "🔻 SHOW " 
	sample_size_2 = parse(Int64, sample_size_2_string)
	E_fraction = expectation_fraction_collected(n, sample_size_2; p_vec = p_vec, r = r)
	
	md""" 	            ↳ `Expected fraction observed`	= **$E_fraction**
	"""	
	end
end

# ╔═╡ 6f839ee1-e980-48f1-a215-f2d040b51823


# ╔═╡ 0099145a-5460-4549-9513-054bc1b04eea
if  show_satur == "🔻 SHOW " 
md""" *A curve describing the expected fraction of modules observed in function of sample size.* """
	end

# ╔═╡ 7968de5e-5ae8-4ab4-b089-c3d33475af2f
begin
	if show_satur == "🔻 SHOW " 
global sample_size_initial_frac = 5
		while (1 - expectation_fraction_collected(n, sample_size_initial_frac; p_vec = p_vec, r = r)) > 0.0005
		global	 sample_size_initial_frac += n/10
		end
	
	sample_sizes_frac = 0: n/10 : sample_size_initial_frac
	
	fracs = expectation_fraction_collected.(n, sample_sizes_frac; p_vec = p_vec, r = r)
	
	plot(sample_sizes_frac, fracs, title = "Expected fraction of modules observed", 
	    xlabel = "sample size", seriestype=:scatter, 
	    ylabel= "E[fraction observed]", label = "", size=(700,400))
end
end

# ╔═╡ 84a2a4de-0667-4120-919a-06e2119112c0


# ╔═╡ f92a6b6e-a556-45cb-a1ae-9f5fe791ffd2
md""" **💻 Occurrence of a specific module**                                                                                                       $(@bind show_occ Select(["🔻 SHOW ", "🔺 HIDE "], default="🔺 HIDE "))\
*How many times one can expect to have collected a specific module in a sample of a given size.*"""

# ╔═╡ ec2a065f-0dc7-44d4-a18b-6c6a228b3ffc
if show_occ == "🔻 SHOW " && distribution != "Zipf's law"
	md"""    👉 Enter the probability of the module of interest: $(@bind p_string TextField(default="0.005"))\
	    👉 Enter the sample size of interest:                $(@bind sample_size_3_string TextField(default="500"))
	""" 	
	
end

# ipv probabiliteit --> rank i: sorteer modules

# ╔═╡ 0e39a993-bb2f-4897-bfe2-5128ec62bef9
if show_occ == "🔻 SHOW " && distribution == "Zipf's law"
	md"""    👉 Enter the rank of the module of interest:        $(@bind rank_string TextField(default="5"))\
	    👉 Enter the sample size of interest:                $(@bind sample_size_4_string TextField(default="500"))
	""" 	
	
end

# ipv probabiliteit --> rank i: sorteer modules

# ╔═╡ 6acb0a97-6469-499f-a5cf-6335d6aa909a
begin

	
if show_occ == "🔻 SHOW " 
	if distribution != "Zipf's law"
	p = parse(Float64, p_string)
	sample_size_3 = parse(Int64, sample_size_3_string)
 	# module_ = 1
# 	p = p_vec[module_]
# 	p = maximum(p_vec) 
	ed = Int(floor(sample_size_3*p))
	j = 0:1:minimum([20, 2*ed])
			
	x  = prob_occurrence_module.(p, sample_size_3, j)
	 plot(j,x, seriestype=[:line, :scatter], xlabel="№ occurrences in sample", ylabel="probability p", title="Chance on № of occurrences for specific module")
	
		else
		rank = parse(Int64, rank_string)
		p = p_vec[rank]
	sample_size_4 = parse(Int64, sample_size_4_string)
 	# module_ = 1
# 	p = p_vec[module_]
# 	p = maximum(p_vec) 
			ed = Int(floor(sample_size_4*p))
	j = 0:1:minimum([20, 2*ed])
			
	x  = prob_occurrence_module.(p, sample_size_4, j)
	 plot(j,x, seriestype=[:line, :scatter], xlabel="№ occurrences in sample", ylabel="probability p", title="Chance on № of occurrences for specific module", size=((550,300)))	
			
		end
	end
end

# ╔═╡ 595423df-728b-43b1-ade4-176785c54be3
begin
	if show_occ == "🔻 SHOW " 
		
	md""" 	            ↳ `Expected number of times observed `	≈ **$ed**
		"""
	end
end

# ╔═╡ fbffaab6-3154-49df-a226-d5810d0b7c38
md"""## References"""

# ╔═╡ 1f48143a-2152-4bb9-a765-a25e70c281a3
md"""[^1]:  Doumas, A. V., & Papanicolaou, V. G. (2016). *The coupon collector’s problem revisited: generalizing the double Dixie cup problem of Newman and Shepp.* ESAIM: Probability and Statistics, 20, 367-399.

[^2]: Boneh, A., & Hofri, M. (1997). *The coupon-collector problem revisited—a survey of engineering problems and computational methods.* Stochastic Models, 13(1), 39-66.



"""


# ╔═╡ Cell order:
# ╟─4d246460-af05-11eb-382b-590e60ba61f5
# ╟─a2fd6000-1450-4dfe-9426-5303ae64bfb3
# ╠═2d3ad982-ef1f-45ae-b247-9679c0faa853
# ╟─a8c81622-194a-443a-891b-bfbabffccff1
# ╟─45507d48-d75d-41c9-a018-299e209f900e
# ╟─b17f3b8a-61ee-4563-97cd-19ff049a8e1e
# ╟─e3b4c2d8-b78c-467e-a863-5eecb8ec58dc
# ╟─2639e3fb-ccbb-44de-bd15-1c5dbf6c1539
# ╟─44d4dfee-3073-49aa-867c-3abea10e6e37
# ╟─f6ebf9fb-0a29-4cb4-a544-6c6e32bedcc4
# ╟─87c3f5cd-79bf-4ad8-b7f8-3e98ec548a9f
# ╟─d877bd4c-497d-46d1-9c58-b6fe26933bfc
# ╟─d4a9da7a-f455-426b-aecd-227c25e1d4e8
# ╟─f098570d-799b-47e2-b692-476a4d95825b
# ╟─caf67b2f-cc2f-4d0d-b619-6e1969fabc1a
# ╟─6f14a72c-51d3-4759-bb8b-10db1dc260f0
# ╟─3c07dd11-5be7-4ec7-992e-38dd07260d71
# ╟─22fe8006-0e81-4e0a-a460-28610a55cd97
# ╟─db4371e4-7f86-4db3-b076-12f6cd220b89
# ╟─317995ed-bdf4-4f78-bd66-a39ffd1dc452
# ╟─9c36cfea-91dd-4bf7-b4f1-b34bce20cb9c
# ╟─ca5a4cef-df67-4a5e-8a86-75a9fe8c6f37
# ╟─24f7aae7-d37a-4db5-ace0-c910b178da88
# ╟─37f951ee-885c-4bbe-a05f-7c5e48ff4b6b
# ╟─ca33610c-1be9-4c01-b0b7-ce4b2f7896df
# ╟─dc696281-7a5b-4568-a4c2-8dde90af43f0
# ╟─eb92ff7c-0140-468c-8b32-f15d1cf15913
# ╟─f0eaf96b-0bc0-4194-9a36-886cb1d66e00
# ╟─6f839ee1-e980-48f1-a215-f2d040b51823
# ╟─0099145a-5460-4549-9513-054bc1b04eea
# ╟─7968de5e-5ae8-4ab4-b089-c3d33475af2f
# ╟─84a2a4de-0667-4120-919a-06e2119112c0
# ╟─f92a6b6e-a556-45cb-a1ae-9f5fe791ffd2
# ╟─ec2a065f-0dc7-44d4-a18b-6c6a228b3ffc
# ╟─0e39a993-bb2f-4897-bfe2-5128ec62bef9
# ╟─6acb0a97-6469-499f-a5cf-6335d6aa909a
# ╟─595423df-728b-43b1-ade4-176785c54be3
# ╟─fbffaab6-3154-49df-a226-d5810d0b7c38
# ╟─1f48143a-2152-4bb9-a765-a25e70c281a3
