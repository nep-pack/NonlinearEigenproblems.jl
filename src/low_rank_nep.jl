export LowRankFactorizedNEP

"
    struct LowRankFactorizedNEP <: AbstractSPMF
    function LowRankFactorizedNEP(L::Vector,U::Vector,f::Vector)
    function LowRankFactorizedNEP(L::Vector,U::Vector,A::Vector, f::Vector)

Representation of a `NEP` which has low rank in the sense that it is an `SPMF`
where each of the terms are factorized: `A[i]=L[i]*U[i]'`. The factorization
is provided in the `L` and `U` vectors and the full matrix `A[i]` can be
either provided (or is otherwise implicitly computed).

# Example:
```
julia> L=randn(5,1); U=randn(5,1);
julia> f=S->exp(S)
julia> nep=LowRankFactorizedNEP([L],[U],[f]);
julia> X=randn(5,2);
julia> norm(compute_Mlincomb(nep,0.0,X)-L*U'*X*ones(2),1)
6.661338147750939e-16
```

"
struct LowRankFactorizedNEP{S<:AbstractMatrix{<:Number}} <: AbstractSPMF{AbstractMatrix}
    spmf::SPMF_NEP
    r::Int          # Sum of ranks of matrices
    L::Vector{S}    # Low rank L factors of matrices
    U::Vector{S}    # Low rank U factors of matrices
end


function LowRankFactorizedNEP(L::AbstractVector{S}, U::AbstractVector{S},
                              A::AbstractVector{S}, f) where {S<:AbstractMatrix}
    rank = mapreduce(u -> size(u, 2), +, U)
    spmf=SPMF_NEP(A, f, align_sparsity_patterns=true, check_consistency=false);
    return LowRankFactorizedNEP{S}(spmf, rank, L, U)
end
function LowRankFactorizedNEP(L::AbstractVector{S}, U::AbstractVector{S},
                              f) where {S<:AbstractMatrix}
    A = L .* adjoint.(U);
    return LowRankFactorizedNEP(L,U,A,f)
end


LowRankFactorizedNEP(::Type{T}, n) where T<:Number =
    LowRankFactorizedNEP(SPMF_NEP(n), 0, Vector{Matrix{T}}(), Vector{Matrix{T}}())

# forward function calls to SPMF
compute_Mder(nep::LowRankFactorizedNEP, λ::T, i::Int = 0) where T<:Number =
    compute_Mder(nep.spmf, λ, i)

compute_Mlincomb(nep::LowRankFactorizedNEP, λ::Number, V::AbstractVecOrMat, a::Vector = ones(eltype(V), size(V, 2)))  =
    compute_Mlincomb(nep.spmf, λ, V, a)

compute_MM(nep::LowRankFactorizedNEP, par...)  =
    compute_MM(nep.spmf, par...)

size(nep::LowRankFactorizedNEP) = size(nep.spmf)
size(nep::LowRankFactorizedNEP, dim) = size(nep.spmf, dim)

get_Av(nep::LowRankFactorizedNEP) = get_Av(nep.spmf)
get_fv(nep::LowRankFactorizedNEP) = get_fv(nep.spmf)
