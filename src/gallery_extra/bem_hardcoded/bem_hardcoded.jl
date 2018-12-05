include("triangle.jl");
include("genmesh.jl");
include("assemble_BEM.jl");
import .NEPCore.compute_Mder
import .NEPCore.size
import .NEPCore.compute_Mlincomb

using LinearAlgebra
"""
    struct BEM_NEP <: NEP

Represents a boundary element discretization. Create with a `mesh::Vector{Triangle}` using `nep=BEM_NEP(mesh)`.

"""
struct BEM_NEP <: NEP
    mesh::Vector{Triangle}
    n::Int
end

function BEM_NEP(mesh::Vector)
    return BEM_NEP(mesh,size(mesh,1));
end

function size(nep::BEM_NEP)
    return (nep.n,nep.n);
end
function size(nep::BEM_NEP,dim)
    return nep.n;
end
#
function compute_Mder(nep::BEM_NEP,λ::Number,der::Int=0)
    return assemble_BEM(λ, nep.mesh, 3,der)[:,:,1];
end
# Delegate the compute Mlincomb functions. This is quite inefficient.
compute_Mlincomb(nep::BEM_NEP,λ::Number,V::Union{AbstractMatrix,AbstractVector}, a::Vector) = compute_Mlincomb_from_Mder(nep,λ,V,a)
compute_Mlincomb(nep::BEM_NEP,λ::Number,V::Union{AbstractMatrix,AbstractVector}) = compute_Mlincomb(nep,λ,V, ones(eltype(V),size(V,2)))



function bem_fichera(N::Int=3)
    mesh=gen_ficheramesh(N)
    return BEM_NEP(mesh);
end
