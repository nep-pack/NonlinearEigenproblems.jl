using NonlinearEigenproblems, Random, SparseArrays, Revise, LinearAlgebra, PyCall, PyPlot
import ..NEPTypes.AbstractSPMF
struct DerSPMF{T<:AbstractMatrix,FDtype,TT<:Number} <: AbstractSPMF{T}
    spmf::SPMF_NEP{T}
    fD::Matrix{FDtype}
    σ::TT
end

# implement all the functions
import ..NEPTypes.size
function size(nep::DerSPMF)
    return size(nep.spmf)
end
function size(nep::DerSPMF,dim)
    return size(nep.spmf,dim)
end

import ..NEPTypes.get_fv
function get_fv(nep::DerSPMF)
    return get_fv(nep.spmf)
end

import ..NEPTypes.get_Av
function get_Av(nep::DerSPMF)
    return get_Av(nep.spmf)
end

## one constructor takes spmf as input and compute the derivatives
function DerSPMF(spmf::SPMF_NEP,σ::Number,m::Int)
      # Compute DD-matrix from get_fv(spmf)
      TT=promote_type(typeof(σ),eltype(nep.A[1]))
      p=length(nep.fi)
      # matrix for the computation of derivatives
      SS=diagm(0=> σ*ones(TT,2m+2),  -1 => (1:2m+1))
      fD=Matrix{TT}(undef, 2*m+2,p)
      for t=1:p fD[:,t]=nep.fi[t](SS)[:,1] end
      return DerSPMF(spmf,fD,σ);
end

import ..NEPCore.compute_Mlincomb
function compute_Mlincomb(
                    nep::DerSPMF{T,FDtype},
                    λ::Number,
                    V::AbstractVecOrMat,
                    a::Vector=ones(size(V,2))) where {T,FDtype}

    local n,k,p
    p=size(nep.fD,2)
    n,k=size(V)
    # Type logic
    TT=promote_type(eltype(V),typeof(λ),eltype(nep.spmf.A[1]),eltype(a))
    z=zeros(TT,n)
    VafD=V*(a.*view(nep.fD,1:k,:));
    @inbounds for j=1:p
        z .+= nep.spmf.A[j]*(view(VafD,:,j))
    end
    return z
end


n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
A4 = sparse(K, J, rand(3*n-2)); A4 = A4+A4';

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
f4= S -> exp(-S)

nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])
σ=0
Dnep=DerSPMF(nep,σ,200)

V,H,ω,HH=ilan(Dnep;Neig=10,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1)
Q=V;

# project (hardcoded for now)
AA1=Q'*(A1*Q);
AA2=Q'*(A2*Q);
AA3=Q'*(A3*Q);
AA4=Q'*(A4*Q);
#err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;
err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;

pnep=SPMF_NEP([AA1,AA2,AA3,AA4],[f1,f2,f3,f4])
λ,_,err=iar(pnep;Neig=100,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);

# sort error
for j=1:p
    err[1:m,j]=sort(err[1:m,j];rev=true);
end

for j=1:p
    semilogy(1:m,err[1:m,j],color="black",linestyle="-");
end
ylim(ymax=1)
