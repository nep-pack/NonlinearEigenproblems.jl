module NleigsTypes

using NEPCore
using NEPTypes

export NleigsNEP
export MatrixAndFunction
export LowRankMatrixAndFunction
export LowRankFactorizedNEP
export NleigsSolutionDetails

import Base.size
import NEPCore.compute_Mder
import NEPCore.compute_Mlincomb
import NEPTypes.get_Av
import NEPTypes.get_fv

struct NleigsNEP{S<:AbstractMatrix{<:Number}, T<:Number}
    nep::NEP            # Original NEP problem
    spmf::Bool          # Whether this NEP is a sum a products of matrices and functions
    p::Int              # Order of polynomial part
    q::Int              # Number of nonlinear matrices and functions
    BBCC::S             # The polynomial and nonlinear matrices concatenated vertically
    is_low_rank::Bool   # Whether the nonlinear matrices of this NEP have an attached low rank LU factorization
    r::Int              # Sum of ranks of low rank nonlinear matrices
    iL::Vector{Int}     # Vector with indices of L-factors
    iLr::Vector{Int}    # Vector with row indices of L-factors that have any entries present
    L::Vector{S}        # Vector of L factors of low rank LU factorization of nonlinear part
    LL::Vector{SparseVector{T,Int}} # Rows of L factors concatenated into vectors (only rows with nnz > 0 are stored)
    UU::S               # The U factors concatenated vertically
end

NleigsNEP(::Type{T}, nep::NEP) where T<:Number =
    NleigsNEP(nep, false, 0, 0, Matrix{T}(0, 0), false, 0, Vector{Int}(0), Vector{Int}(0), Vector{Matrix{T}}(0), Vector{SparseVector{T,Int}}(0), Matrix{T}(0, 0))

NleigsNEP(nep::NEP, p, q, BBCC::AbstractMatrix{T}) where T<:Number =
    NleigsNEP(nep, true, p, q, BBCC, false, 0, Vector{Int}(0), Vector{Int}(0), Vector{Matrix{T}}(0), Vector{SparseVector{T,Int}}(0), Matrix{T}(0, 0))

NleigsNEP(nep::NEP, p, q, BBCC, r, iL, iLr, L, LL, UU) =
    NleigsNEP(nep, true, p, q, BBCC, true, r, iL, iLr, L, LL, UU)

struct LowRankMatrixAndFunction{S<:AbstractMatrix{<:Number}}
    A::S
    L::S        # L factor of LU-factorized A
    U::S        # U factor of LU-factorized A
    f::Function
end

"Create low rank LU factorization of A."
function LowRankMatrixAndFunction(A::AbstractMatrix{<:Number}, f::Function)
    L, U = low_rank_lu_factors(A)
    LowRankMatrixAndFunction(A, L, U, f)
end

function low_rank_lu_factors(A::SparseMatrixCSC{<:Number,Int64})
    n = size(A, 1)
    r, c = findn(A)
    r = extrema(r)
    c = extrema(c)
    B = A[r[1]:r[2], c[1]:c[2]]
    L, U = lu(full(B))
    Lc, Uc = compactlu(sparse(L), sparse(U))
    Lca = spzeros(n, size(Lc, 2))
    Lca[r[1]:r[2], :] = Lc
    Uca = spzeros(size(Uc, 1), n)
    Uca[:, c[1]:c[2]] = Uc
    return Lca, Uca'

    # TODO use this; however we then need to support permutation and scaling
    #F = lufact(B)
    #Lcf,Ucf = compactlu(sparse(F[:L]),sparse(F[:U]))
    #Lcaf = spzeros(n, size(Lcf, 2))
    #Lcaf[r[1]:r[2], :] = Lcf
    #Ucaf = spzeros(size(Ucf, 1), n)
    #Ucaf[:, c[1]:c[2]] = Ucf
    #Ucaf = Ucaf'
end

function compactlu(L, U)
    n = size(L, 1)
    select = map(i -> nnz(L[i:n, i]) > 1 || nnz(U[i, i:n]) > 0, 1:n)
    return L[:,select], U[select,:]
end

"SPMF with low rank LU factors for each matrix."
struct LowRankFactorizedNEP{S<:AbstractMatrix{<:Number}} <: AbstractSPMF
    spmf::SPMF_NEP
    r::Int          # Sum of ranks of matrices
    L::Vector{S}    # Low rank L factors of matrices
    U::Vector{S}    # Low rank U factors of matrices
end

function LowRankFactorizedNEP(Amf::AbstractVector{LowRankMatrixAndFunction{S}}) where {T<:Number, S<:AbstractMatrix{T}}
    q = length(Amf)
    r = 0
    f = Vector{Function}(q)
    A = Vector{S}(q)
    L = Vector{S}(q)
    U = Vector{S}(q)

    for k = 1:q
        f[k] = Amf[k].f
        L[k] = Amf[k].L
        U[k] = Amf[k].U
        # if A is not specified, create it from LU factors
        A[k] = isempty(Amf[k].A) ? L[k] * U[k]' : Amf[k].A
        r += size(U[k], 2)
    end

    return LowRankFactorizedNEP(SPMF_NEP(A, f), r, L, U)
end

"Create an empty LowRankFactorizedNEP."
LowRankFactorizedNEP(::Type{T}, n) where T<:Number =
    LowRankFactorizedNEP(SPMF_NEP(n), 0, Vector{Matrix{T}}(0), Vector{Matrix{T}}(0))

# forward function calls to SPMF
compute_Mder(nep::LowRankFactorizedNEP, λ::T, i::Int = 0) where T<:Number =
    compute_Mder(nep.spmf, λ, i)

compute_Mlincomb(nep::LowRankFactorizedNEP, λ::T, V::Union{Vector{T}, Matrix{T}}; a = ones(T, size(V, 2))) where T<:Number =
    compute_Mlincomb(nep.spmf, λ, V, a=a)

size(nep::LowRankFactorizedNEP) = size(nep.spmf)
size(nep::LowRankFactorizedNEP, dim) = size(nep.spmf, dim)

get_Av(nep::LowRankFactorizedNEP) = get_Av(nep.spmf)
get_fv(nep::LowRankFactorizedNEP) = get_fv(nep.spmf)

struct NleigsSolutionDetails{T<:Real, CT<:Complex{T}}
    "matrix of Ritz values in each iteration"
    Lam::AbstractMatrix{CT}

    "matrix of residuals in each iteraion"
    Res::AbstractMatrix{T}

    "vector of interpolation nodes"
    sigma::AbstractVector{CT}

    "vector of poles"
    xi::AbstractVector{T}

    "vector of scaling parameters"
    beta::AbstractVector{T}

    "vector of norms of generalized divided differences (in function handle
    case) or maximum of absolute values of scalar divided differences in
    each iteration (in matrix function case)"
    nrmD::AbstractVector{T}

    "number of iterations until linearization converged"
    kconv::Int
end

NleigsSolutionDetails{T,CT}() where {T<:Real, CT<:Complex{T}} = NleigsSolutionDetails(
    Matrix{CT}(0,0), Matrix{T}(0, 0), Vector{CT}(0),
    Vector{T}(0), Vector{T}(0), Vector{T}(0), 0)

end
