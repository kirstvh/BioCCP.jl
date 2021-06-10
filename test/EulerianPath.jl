@testset "EulerianPath" begin
    using STMOZOO.EulerianPath
    adj_list = create_adj_list([[1, 5], [1, 2], [1, 4], [1, 3], [5, 2], [2, 3], [2, 4], [3, 4], [3, 6], [4, 6]])
    @testset "create_adj_list" begin
        
        @test adj_list isa Dict{Int64,Array{Int64,N} where N}
        @test adj_list == Dict(4 => [1, 2, 3, 6],2 => [1, 5, 3, 4],3 => [1, 2, 4, 6],5 => [1, 2],6 => [3, 4],1 => [5, 2, 4, 3])

    end

    @testset "has_eulerian_cycle" begin
        @test has_eulerian_cycle(adj_list) == true
    end

    @testset "has_eulerian_path" begin
        @test has_eulerian_path(adj_list) == true
    end

    adj_list = create_adj_list([[1, 2], [1, 4], [1, 3], [2, 3], [2, 4], [3, 4]])
    @testset "create_adj_list" begin
        @test adj_list == Dict(4 => [1, 2, 3],2 => [1, 3, 4],3 => [1, 2, 4],1 => [2, 4, 3])
    end

    @testset "has_eulerian_cycle" begin
        @test has_eulerian_cycle(adj_list) == false
    end

    @testset "has_eulerian_path" begin
        @test has_eulerian_path(adj_list) == false
    end

end