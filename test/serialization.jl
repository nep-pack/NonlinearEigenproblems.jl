workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src/utils"))

using Base.Test
using Serialization

@testset "Serialization" begin
    file = "sparse_matrix.txt"

    A = sprand(40, 20, 0.05)
    write_sparse_matrix(file, A)

    B = read_sparse_matrix(file)

    rm(file)

    @test A == B
end
