export ConvergedLinearization
export Linearization
export CorkSolutionDetails

struct ConvergedLinearization
    n::Int # Size of problem
    p::Int # Order of polynomial part
    U::AbstractMatrix # U factors of the low rank nonlinear part
    d::Int # Degree of approximation
    σ::Vector # Interpolation nodes
    ξ::Vector # Poles
    β::Vector # Scaling factors
    D::Vector #{<:AbstractMatrix} # Full rank generalized divided differences
    Dt::Vector # Low rank generalized divided differences
end

"""
A:     cell array of n x n matrices of length d
B:     cell array of n x n matrices of length d
M:     (d-1) x d matrix
N:     (d-1) x d matrix

      [ A_0  A_1  ...  A_{d-1} ]         [ B_0  B_1  ...  B_{d-1} ]
      [ ---------------------- ]         [ ---------------------- ]
L  =  [                        ] - lam * [                        ]
      [     kron(M,eye(n))     ]         [     kron(N,eye(n))     ]
      [                        ]         [                        ]
"""
struct Linearization
    A::Vector # vector of n x n matrices of length d
    B::Vector # vector of n x n matrices of length d
    M::AbstractMatrix # (d-1) x d matrix
    N::AbstractMatrix # (d-1) x d matrix
end

struct CorkSolutionDetails{T<:Real, CT<:Complex{T}}
    "matrix of Ritz values in each iteration"
    Lam::Matrix{CT}

    "matrix of residuals in each iteraion"
    Res::Matrix{T}

    "vector with dimension of Krylov subspace in each iteration"
    J::Vector{Int}

    "vector with rank of subspace Q in each iteration"
    R::Vector{Int}

    "maximum number of ritz values"
    m::Int

    "number of restarted ritz values"
    p::Int
end

CorkSolutionDetails{T,CT}() where {T<:Real, CT<:Complex{T}} = CorkSolutionDetails(
    Matrix{CT}(undef, 0, 0), Matrix{T}(undef, 0, 0), Vector{Int}(), Vector{Int}(), 0, 0)

include("method_cork.jl")
