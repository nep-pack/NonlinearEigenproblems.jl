using NonlinearEigenproblems
using SparseArrays

export RKNEP
export MatrixAndFunction
export LowRankMatrixAndFunction
export LowRankFactorizedNEP

export get_rk_nep

import Base.size
import ..NEPCore.compute_Mder
import ..NEPCore.compute_Mlincomb
import ..NEPTypes.get_Av
import ..NEPTypes.get_fv

"NEP instance supplemented with various structures used for rational Krylov problems."
struct RKNEP{S<:AbstractMatrix{<:Number}, T<:Number}
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

RKNEP(::Type{T}, nep::NEP) where T<:Number =
    RKNEP(nep, false, 0, 0, Matrix{T}(undef, 0, 0), false, 0, Vector{Int}(), Vector{Int}(), Vector{Matrix{T}}(), Vector{SparseVector{T,Int}}(), Matrix{T}(undef, 0, 0))

RKNEP(nep::NEP, p, q, BBCC::AbstractMatrix{T}) where T<:Number =
    RKNEP(nep, true, p, q, BBCC, false, 0, Vector{Int}(), Vector{Int}(), Vector{Matrix{T}}(), Vector{SparseVector{T,Int}}(), Matrix{T}(undef, 0, 0))

RKNEP(nep::NEP, p, q, BBCC, r, iL, iLr, L, LL, UU) =
    RKNEP(nep, true, p, q, BBCC, true, r, iL, iLr, L, LL, UU)

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
    idx = findall(!iszero, A)
    r = extrema(getindex.(idx, 1))
    c = extrema(getindex.(idx, 2))
    B = A[r[1]:r[2], c[1]:c[2]]
    L, U = lu(Matrix(B); check = false)
    Lc, Uc = compactlu(sparse(L), sparse(U))
    Lca = spzeros(n, size(Lc, 2))
    Lca[r[1]:r[2], :] = Lc
    Uca = spzeros(size(Uc, 1), n)
    Uca[:, c[1]:c[2]] = Uc
    return Lca, sparse(Uca')

    # TODO use this; however we then need to support permutation and scaling
    #F = lu(B)
    #Lcf,Ucf = compactlu(sparse(F.L),sparse(F.U))
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
struct LowRankFactorizedNEP{S<:AbstractMatrix{<:Number}} <: AbstractSPMF{AbstractMatrix}
    spmf::SPMF_NEP
    r::Int          # Sum of ranks of matrices
    L::Vector{S}    # Low rank L factors of matrices
    U::Vector{S}    # Low rank U factors of matrices
end

function LowRankFactorizedNEP(Amf::AbstractVector{LowRankMatrixAndFunction{S}}) where {T<:Number, S<:AbstractMatrix{T}}
    q = length(Amf)
    r = 0
    f = Vector{Function}(undef, q)
    A = Vector{S}(undef, q)
    L = Vector{S}(undef, q)
    U = Vector{S}(undef, q)

    for k = 1:q
        f[k] = Amf[k].f
        L[k] = Amf[k].L
        U[k] = Amf[k].U
        # if A is not specified, create it from LU factors
        A[k] = isempty(Amf[k].A) ? L[k] * U[k]' : Amf[k].A
        r += size(U[k], 2)
    end

    return LowRankFactorizedNEP(SPMF_NEP(A, f, align_sparsity_patterns=true), r, L, U)
end

"Create an empty LowRankFactorizedNEP."
LowRankFactorizedNEP(::Type{T}, n) where T<:Number =
    LowRankFactorizedNEP(SPMF_NEP(n), 0, Vector{Matrix{T}}(), Vector{Matrix{T}}())

# forward function calls to SPMF
compute_Mder(nep::LowRankFactorizedNEP, λ::T, i::Int = 0) where T<:Number =
    compute_Mder(nep.spmf, λ, i)

compute_Mlincomb(nep::LowRankFactorizedNEP, λ::T, V::Union{Vector{T}, Matrix{T}}, a::Vector = ones(T, size(V, 2))) where T<:Number =
    compute_Mlincomb(nep.spmf, λ, V, a)

size(nep::LowRankFactorizedNEP) = size(nep.spmf)
size(nep::LowRankFactorizedNEP, dim) = size(nep.spmf, dim)

get_Av(nep::LowRankFactorizedNEP) = get_Av(nep.spmf)
get_fv(nep::LowRankFactorizedNEP) = get_fv(nep.spmf)

"Create RKNEP instance, exploiting the type of the input NEP as much as possible"
function get_rk_nep(::Type{T}, nep::NEP) where T<:Real
    # Most generic case: No coefficient matrices, all we have is M(λ)
    if !isa(nep, AbstractSPMF)
        return RKNEP(T, nep)
    end

    Av = get_Av(nep)
    BBCC = vcat(Av...)::eltype(Av)

    # Polynomial eigenvalue problem
    if isa(nep, PEP)
        return RKNEP(nep, length(Av) - 1, 0, BBCC)
    end

    # If we can't separate the problem into a PEP + SPMF, consider it purely SPMF
    if !isa(nep, SPMFSumNEP{PEP,S} where S<:AbstractSPMF)
        return RKNEP(nep, -1, length(Av), BBCC)
    end

    p = length(get_Av(nep.nep1)) - 1
    q = length(get_Av(nep.nep2))

    # Case when there is no low rank structure to exploit
    if q == 0 || !isa(nep.nep2, LowRankFactorizedNEP{S} where S<:Any)
        return RKNEP(nep, p, q, BBCC)
    end

    # L and U factors of the low rank nonlinear part
    L = nep.nep2.L
    UU = hcat(nep.nep2.U...)::eltype(nep.nep2.U)
    r = nep.nep2.r
    iL = zeros(Int, r)
    c = 0
    for ii = 1:q
        ri = size(L[ii], 2)
        iL[c+1:c+ri] .= ii
        c += ri
    end

    # Store L factors in a compact format to speed up system solves later on
    LL = Vector{SparseVector{eltype(L[1]),Int}}()
    iLr = Vector{Int}()
    for ri = 1:size(nep, 1)
        row = reduce(vcat, [L[i][ri,:] for i=1:length(L)])
        if nnz(row) > 0
            push!(LL, row)
            push!(iLr, ri)
        end
    end

    return RKNEP(nep, p, q, BBCC, r, iL, iLr, L, LL, UU)
end
