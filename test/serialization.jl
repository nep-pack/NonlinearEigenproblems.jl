push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems.Serialization
using Test
using SparseArrays

@bench @testset "Serialization" begin
    file = "sparse_matrix.txt"

    A = sprand(40, 20, 0.05)
    write_sparse_matrix(file, A)

    B = read_sparse_matrix(file)

    rm(file)

    @test A == B
end
