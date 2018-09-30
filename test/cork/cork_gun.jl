using NonlinearEigenproblems

include("../nleigs/gun_test_utils.jl")
include("../../src/cork/linearize.jl")

#function cork_gun()
    nep, Σ, Ξ, _, nodes, funres = gun_init()

    CL = linearize(Float64, nep, Σ, Ξ, 100, 1e-10)

    n = CL.n
    d = CL.d

    A = CL.D
    B = Vector{SparseMatrixCSC}(undef, d)
    for i = 1:d
        B[i] = spzeros(n, n)
    end
    CT = eltype(CL.σ)
    M = [diagm(0 => CL.σ[1:d-1], 1 => CL.β[2:d-1]) [zeros(CT, d-2); CL.β[d]]]
    N = [I + diagm(1 => CL.β[2:d-1]./CL.ξ[1:d-2]) [zeros(CT, d-2); CL.β[d]/CL.ξ[d-1]]]

    L = Linearization(A, B, M, N)

    Random.seed!(0)
    v = randn(n) .+ im * randn(n) # TODO: use gun_init v?

    @time lambda, X, res, flag, solution_info =
        cork(Float64, nep, L, nodes[3], Σ, k=20, m=50, p=35, tolres=1e-10, shifts=nodes, v0=v, errmeasure=funres, displaylevel=1)

    println("Found $(length(lambda)) lambdas:")
    foreach(x -> println("  $x"), lambda)

#end

#cork_gun()
