using Base: Int64
using BioCCP
@testset "BioCCP" begin
    
    n = 100

    # Equal probabilities
    p_uniform = ones(n)/n

    # Probabilities following Zipf's law
    ρ = 10
	α = exp(log(ρ)/(n-1))
	p_zipf = collect(α.^-(1:n))
    p_zipf = p_zipf ./ sum(p_zipf)

    @testset "Expectation" begin
        @test expectation_minsamplesize(n; p = p_uniform) isa Int64
        @test expectation_minsamplesize(n; p = p_uniform) < expectation_minsamplesize(n; p = p_zipf)
    end

    @testset "Standard deviation" begin
        @test std_minsamplesize(n; p = p_uniform) isa Int64
        @test std_minsamplesize(n; p = p_uniform) < std_minsamplesize(n; p = p_zipf)
    end 

    t = 500

    @testset "Success probability" begin
        @test success_probability(n, t; p = p_uniform) isa Float64
        @test success_probability(n, t; p = p_uniform) > success_probability(n, t; p = p_zipf)
    end

    @testset "Expected fraction collected" begin
        @test expectation_fraction_collected(n, t; p = p_uniform) isa Float64
        @test expectation_fraction_collected(n, t; p = p_uniform) > expectation_fraction_collected(n, t; p = p_zipf)
    end

    k = 2
    pᵢ = 0.005

    @testset "Probability occurence" begin
        @test prob_occurrence_module(pi, t, k) isa Float64
    end
end