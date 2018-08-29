if !@isdefined global_modules_loaded
    using Test
end

include(normpath(string(@__DIR__), "..", "src", "utils", "Serialization.jl"))

@testset "Serialization" begin
    file = "sparse_matrix.txt"

    A = sprand(40, 20, 0.05)
    write_sparse_matrix(file, A)

    B = read_sparse_matrix(file)

    rm(file)

    @test A == B
end
