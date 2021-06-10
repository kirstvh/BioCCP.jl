@testset "MaximumFlow" begin
    using STMOZOO.MaximumFlow,LinearAlgebra
    test_mat = [0 8 2 0; 0 0 2 4;0 0 0 1;0 0 0 0] # case
    test_mat2 = [0 7 0 3 0; 0 0 2 0 2; 0 0 0 0 4; 0 0 1 0 2; 0 0 0 0 0] # case
    cur_mat = [0 1 1 0; 0 0 1 1;0 0 0 1;0 0 0 0] # case
    test_mult = [0 0 0 4 0 0 0 0; 0 0 0 7 0 0 0 0; 0 0 0 0 2 0 0 0; 0 0 0 0 5 0 5 10;
                 0 0 0 0 0 8 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0] # case

    @testset "res_network" begin
        @test res_network(test_mat, cur_mat) isa Array{Int,2} # test type
        @test res_network(test_mat, test_mat) == test_mat' # test case
        @test all(res_network(test_mat, cur_mat) .>= 0) # # test property
        @test tr(res_network(test_mat, cur_mat)) == 0 # # test property
        @test res_network(test_mat, zeros(Int, 4, 4)) == test_mat # test case
        @test_throws AssertionError res_network(rand(Int, 3, 3), rand(Int, 3, 4)) # test error
        @test_throws AssertionError res_network(rand(Int, 5, 3), rand(Int, 4, 4)) # test error
        @test_throws AssertionError res_network(rand(Int, 5, 5), rand(Int, 4, 4)) # test error
    end

    @testset "augmenting_path" begin
        flow_a, mflow_a = augmenting_path(test_mat, 1, size(test_mat,1)) # case
        flow_b, mflow_b = augmenting_path(test_mat, 1, 2) # case
        flow_c, mflow_c = augmenting_path(test_mat, 1, 3) # case
        flow2_a, mflow2_a = augmenting_path(test_mat2, 1, size(test_mat2,2)) # case
        flow2_b, mflow2_b = augmenting_path(test_mat2, 1, 2) # case
        flow2_c, mflow2_c = augmenting_path(test_mat2, 1, 3) # case
        flow2_d, mflow2_d = augmenting_path(test_mat2, 1, 4) # case
        @test flow_a isa Array{Int,2} # test type
        @test mflow_a isa Int # test type
        @test tr(flow_a) == 0 # # test property
        @test all([sum(flow_a[i,:]) == sum(flow_a[:,i]) for i in 2:(size(test_mat,1)-1)]) # test property
        @test all([sum(flow2_a[i,:]) == sum(flow2_a[:,i]) for i in 2:(size(test_mat2,1)-1)]) # test property
        @test sum(flow_a[:,size(test_mat,1)]) == sum(flow_a[1,:]) # test property
        @test mflow_a == 5 # test case
        @test mflow_b == 8 # test case
        @test mflow_c == 4 # test case
        @test mflow2_a == 7 # test case
        @test mflow2_b == 7 # test case
        @test mflow2_c == 3 # test case
        @test mflow2_d == 3 # test case
        @test_throws AssertionError augmenting_path(rand(Int, 5, 3), 1,3) # test error
        @test augmenting_path(test_mult,1,8)[2] == 4 # test case
        @test augmenting_path(test_mult,[1, 2, 3],8)[2] == 10 # test case
        @test augmenting_path(test_mult,2,[6, 7, 8])[2] == 7 # test case
        @test augmenting_path(test_mult,[1, 2, 3],[6, 7 ,8])[2] == 13 # test case
        @test augmenting_path(test_mult,[1 ,2],[8 ,7])[2] == 11 # test case

    end
    
    @testset "bfs" begin
        no_path = [0 10 4 0; 7 0 4 0;0 5 0 0;0 0 0 0] # case
        @test bfs(test_mat, 1, 4)[2] isa Array{Int,1} # test type
        @test all(bfs(test_mat, 1, 4)[2] .>= 0) # test property
        @test bfs(test_mat, 1, 4)[1] # test case
        @test bfs(test_mat2, 1, 5)[1] # test case
        @test ! bfs(no_path, 1, 4)[1] # test case
    end

    @testset "isfeasible" begin
        flow_a,mflow_a = augmenting_path(test_mult,[1,2,3],[6,7,8], send = [3,3,2]) # case
        flow_b,mflow_b = augmenting_path(test_mult,[1,2,3],[6,7,8], desired = [1,1,1]); # case
        flow_c,mflow_c = augmenting_path(test_mult,[1,3],[6,7,8], desired = [1,1,1]); # case
        flow_d,mflow_d = augmenting_path(test_mult,[1],[8], send = [3],desired = [3]); # case
        flow_e,mflow_e = augmenting_path(test_mult,[1,2],[8], send = [4,6],desired = [10]); # case
        flow_g,mflow_g = augmenting_path(test_mult,[1,3],[6,8], send = [3,1],desired = [1,3]); # case
        @test isfeasible(flow_a,s=[1,2,3],t=[6,7,8], send = [3,3,2]) # test case
        @test isfeasible(flow_b,s=[1,2,3],t=[6,7,8], desired = [1,1,1]) # test case
        @test isfeasible(flow_c,s=[1,3],t=[6,7,8], desired = [1,1,1]) # test case
        @test isfeasible(flow_d,s=[1],t=[8], send = [3],desired = [3]) # test case
        @test isfeasible(flow_e,s=[1,2],t=[8], send = [4,6],desired = [10]) # test case
        @test isfeasible(flow_g,s=[1,3],t=[6,8], send = [3,1],desired = [1,3]) # test case
        @test_throws AssertionError isfeasible(flow_e,s=[1],t=[8], send = [4,6],desired = [10])
        @test_throws AssertionError isfeasible(flow_e,s=[1,2],t=[7,8], send = [4,6],desired = [10])
    end

    @testset "clean!" begin
        flow_a,mflow_a = augmenting_path(test_mult,[1,2],[8], send = [4,6],desired = [10]); # case
        flow_b,mflow_b = augmenting_path(test_mult,[1,3],[6,8], send = [3,1],desired = [1,3]); # case
        clean!(flow_a)
        clean!(flow_b)
        @test all(flow_a .>=0)
        @test all(flow_b .>=0)
        @test isfeasible(flow_a,s=[1,2],t=[8], send = [4,6],desired = [10]) # test case
        @test isfeasible(flow_b,s=[1,3],t=[6,8], send = [3,1],desired = [1,3]) # test case
    end
end