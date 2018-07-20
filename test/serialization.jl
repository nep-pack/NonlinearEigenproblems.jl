workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src/utils"))

using Base.Test
using Serialization

@testset "Serialization" begin
    file = "sparse_matrix.zip"

    A1 = sprand(40, 20, 0.05)
    A2 = sprand(20, 30, 0.1)
    write_sparse_matrices(file, Dict("A1" => A1, "A2" => A2))

    D = read_sparse_matrices(file)

    rm(file)

    @test A1 == D["A1"]
    @test A2 == D["A2"]
end
