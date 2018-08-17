module NleigsTypes

using NEPTypes

export MatrixAndFunction
export LowRankMatrixAndFunction
export LowRankFactorizedNEP

import Base.size
import NEPCore.compute_Mder
import NEPCore.compute_Mlincomb

abstract type AbstractMatrixAndFunction{S<:AbstractMatrix{<:Real}} end

struct MatrixAndFunction{S<:AbstractMatrix{<:Real}} <: AbstractMatrixAndFunction{S}
    A::S
    f::Function
end

struct LowRankMatrixAndFunction{S<:AbstractMatrix{<:Real}} <: AbstractMatrixAndFunction{S}
    A::S
    L::S        # L factor of LU-factorized A
    U::S        # U factor of LU-factorized A
    f::Function
end

"Creates low rank LU factorization of A"
function LowRankMatrixAndFunction(A::AbstractMatrix{<:Real}, f::Function)
    L, U = low_rank_lu_factors(A)
    LowRankMatrixAndFunction(A, L, U, f)
end

function low_rank_lu_factors(A::SparseMatrixCSC{<:Real,Int64})
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

"""
SPMF with low rank LU factors for each matrix.
"""
struct LowRankFactorizedNEP{S<:AbstractMatrix{<:Real}} <: AbstractSPMF
    spmf::SPMF_NEP
    r::Int          # Sum of ranks of matrices
    L::Vector{S}    # Low rank L factors of matrices
    U::Vector{S}    # Low rank U factors of matrices
end

function LowRankFactorizedNEP(Amf::AbstractVector{<:AbstractMatrixAndFunction{S}}) where {T<:Real, S<:AbstractMatrix{T}}
    q = length(Amf)
    r = 0
    f = Vector{Function}(q)
    A = Vector{S}(q)
    # L and U factors of low rank A
    L = Vector{S}(0)
    U = Vector{S}(0)

    if q > 0
        for k = 1:q
            f[k] = Amf[k].f
            A[k] = Amf[k].A
        end

        if isa(Amf[1], LowRankMatrixAndFunction)
            L = Vector{S}(q)
            U = Vector{S}(q)
            for k = 1:q
                L[k] = Amf[k].L
                U[k] = Amf[k].U
                r += size(U[k], 2)
                # if A is not specified, create it from LU factors
                isempty(A[k]) && (A[k] = L[k] * U[k]')
            end
        end
    end
    return LowRankFactorizedNEP(SPMF_NEP(A, f), r, L, U)
end

# Create an empty LowRankFactorizedNEP
LowRankFactorizedNEP(::Type{T}, n) where T<:Number =
    LowRankFactorizedNEP(SPMF_NEP(n), 0, Vector{Matrix{T}}(0), Vector{Matrix{T}}(0))

# forward function calls to SPMF
compute_Mder(nep::LowRankFactorizedNEP, λ::T, i::Int = 0) where T<:Number =
    compute_Mder(nep.spmf, λ, i)

compute_Mlincomb(nep::LowRankFactorizedNEP, λ::T, V::Union{Vector{T}, Matrix{T}}; a = ones(T, size(V, 2))) where T<:Number =
    compute_Mlincomb(nep.spmf, λ, V, a=a)

size(nep::LowRankFactorizedNEP) = size(nep.spmf)
size(nep::LowRankFactorizedNEP, dim) = size(nep.spmf, dim)

end
