using NonlinearEigenproblems.Serialization
using Test
using SparseArrays

@bench @testset "Serialization" begin
    file = "sparse_matrix.txt"

    Random.seed!(0)
    A = sprand(40, 20, 0.05)
    write_sparse_matrix(file, A)

    B = read_sparse_matrix(file)

    rm(file)

    @test A == B
end
