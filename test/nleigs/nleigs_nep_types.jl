# Solves a basic eigenvalue problem defined through different NEP types through NLEIGS

#using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using NonlinearEigenproblems.RKHelper
using Test
using SparseArrays

struct CustomNLEIGSNEP <: NEP
    n::Int  # problem size; this is the only required field in a custom NEP type when used with NLEIGS
end

nep_types_test_n = 2
nep_types_test_B = Vector{Matrix{Float64}}([[1 3; 5 6], [3 4; 6 6]])
nep_types_test_C = Vector{Matrix{Float64}}([[1 0; 0 1]])

# implement a few methods used by the solver
import NonlinearEigenproblems.NEPCore.compute_Mder
import NonlinearEigenproblems.NEPCore.compute_Mlincomb
import Base.size
nep_types_test_pep = PEP([nep_types_test_B; nep_types_test_C])
compute_Mder(::CustomNLEIGSNEP, λ::Number) = compute_Mder(nep_types_test_pep, λ)
compute_Mlincomb(::CustomNLEIGSNEP, λ::Number, x::AbstractVecOrMat) = compute_Mlincomb(nep_types_test_pep, λ, x)
size(::CustomNLEIGSNEP, _) = nep_types_test_n

function nleigs_nep_types()
    n = nep_types_test_n
    B = nep_types_test_B
    C = nep_types_test_C
    f = [λ -> λ^2]
    Σ = [-10.0-2im, 10-2im, 10+2im, -10+2im]

    # define and solve the same problem in many different ways
    problems = [
        ("SPMF_NEP", SPMF_NEP([B; C], [λ -> λ^0; λ -> λ; λ -> λ^2])),
        ("PEP", PEP([B; C])),
        ("PEP + SPMF", SumNEP(PEP(B), SPMF_NEP(C, f))),
        ("PEP + LowRankFactorizedNEP", SumNEP(PEP(B), LowRankFactorizedNEP([LowRankMatrixAndFunction(sparse(C[1]), f[1])]))),
        ("Custom NEP type", CustomNLEIGSNEP(n))]

    for problem in problems
        @bench @testset "$(problem[1])" begin
            lambda, X = nleigs(problem[2], Σ, maxit=10, v=ones(n).+0im, blksize=5)
            verify_lambdas(4, problem[2], lambda, X)
        end
    end
end

@testset "NLEIGS: NEP Types" begin
    nleigs_nep_types()
end
