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
        @test expectation_minsamplesize(n; p = p_uniform) == Int(ceil(n*sum(1 ./ (1:n))))
        @test expectation_minsamplesize(n; p = p_uniform) < expectation_minsamplesize(n; p = p_zipf)
        @test expectation_minsamplesize(n; p = p_uniform, m = 1) < expectation_minsamplesize(n; p = p_uniform, m = 2)
        @test expectation_minsamplesize(n; p = p_uniform, m = 1, r = 1)/3 == expectation_minsamplesize(n; p = p_uniform, m = 1, r = 3)
    end

    @testset "Standard deviation" begin
        @test std_minsamplesize(n; p = p_uniform) isa Int64
        @test std_minsamplesize(n; p = p_uniform) < std_minsamplesize(n; p = p_zipf)
        @test std_minsamplesize(n; p = p_uniform, m = 1) < std_minsamplesize(n; p = p_uniform, m = 2)
    end 

    t = 500

    @testset "Success probability" begin
        @test success_probability(n, t; p = p_uniform) isa Float64
        @test success_probability(n, t; p = p_uniform) > success_probability(n, t; p = p_zipf)
        @test success_probability(n, t; p = p_uniform, m = 1) > success_probability(n, t; p = p_uniform, m = 2)
        @test all(success_probability(n, t; p = p_uniform) .< success_probability.(n, [t+10, t+100, t+1000]; p = p_uniform))
        @test all(success_probability.(n, [0, 10, 100, 1000, 10000]) .<= 1)
        @test all(success_probability.(n, [0, 10, 100, 1000, 10000]) .>= 0)

    end

    @testset "Expected fraction collected" begin
        @test expectation_fraction_collected(n, t; p = p_uniform) isa Float64
        @test expectation_fraction_collected(n, t; p = p_uniform) > expectation_fraction_collected(n, t; p = p_zipf)
        @test all(expectation_fraction_collected(n, t; p = p_uniform) .< expectation_fraction_collected.(n, [t+10, t+100, t+1000]; p = p_uniform))
        @test all(expectation_fraction_collected.(n, [0, 10, 100, 1000, 10000]) .<= 1)
        @test all(expectation_fraction_collected.(n, [0, 10, 100, 1000, 10000]) .>= 0)

    end

    k = 2
    pᵢ = 0.005
    r = 3

    @testset "Probability occurence" begin
        @test prob_occurrence_module(pᵢ, t, r, k) isa Float64
        @test all(prob_occurrence_module.(pᵢ, [0, 10, 100, 1000, 10000], r, k) .>= 0)
        @test all(prob_occurrence_module.(pᵢ, [0, 10, 100, 1000, 10000], r, k) .<= 1)
    end
end
