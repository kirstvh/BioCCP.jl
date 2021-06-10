@testset "GenProgAlign" begin
using STMOZOO.GenProgAlign
    @testset "ga_alignment" begin
    raw_sequences, seq_ids_raw = fasta_to_array("./notebook/UBQ_small.fasta")
    alignment, score, seq_ids = ga_alignment("./notebook/UBQ_small.fasta", population_size = 10, p = 0.95, max_operations_without_improvement = 50, max_operations = 2000)
        @test isa(alignment, Array{String, 1})
        @test isa(score, Float64)
        @test isa(seq_ids, Array{String, 1})
        @test all([replace(alignment[i], "-" => "") == raw_sequences[i] for i in 1:10])
        @test all(seq_ids_raw .== seq_ids)
        @test all([length(alignment[seq]) for seq in 1:10] .== length(alignment[1]))
        @test  all([alignment[i][j] in vcat(AA, '-') for i = 1:10 for j = 1:length(alignment[1])])
    end
end
