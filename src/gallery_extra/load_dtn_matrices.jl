# Script which loads the DtN nep.
import Base.size
using Printf, SpecialFunctions;
using NonlinearEigenproblems;

include("petsc_naive_bin_read.jl")

l=40; # Number of terms

data_dir="/home/jarl/archive/dtn_umea_collab_nosync_julia/dimer_TM_ref2_p10_g0_a1"

A = naive_petsc_read(joinpath(data_dir,"K.bin"))
M = naive_petsc_read(joinpath(data_dir,"M.bin"))
n = size(A,1);

# config_fem.m contains dofs and dofs_boundaries.
# This should be the same as dof-dofs_boundaries
q1 = naive_petsc_read(joinpath(data_dir,"q1.bin"))
start_dtn_dofs=findfirst(abs.(q1).>0)-1;

dtn_a = 1;

# Determine the number of saved terms.
searchdir(path,key) = filter(x->startswith(x,key), readdir(path))
files=searchdir(data_dir,"q") # All the Q-vectors
# Mid of the number of saved terms.
mid=round(Int,(size(files,1)-1)/2+1);

l=min((mid-1),l); # make sure we don't exceed the dimension

ind=Vector((mid-l):(mid+l));
ind2=ind .- mid;

Q=Matrix{ComplexF64}(undef,n,size(ind,1)); # Factorized matrices
P=Vector{SparseMatrixCSC}(undef,size(ind,1)) # A vector of rank one matrices
for i=1:size(ind,1)
    # This corresponds to the term q*q'*bessel_quotient(ind2[i],s)
    # Load the q-vector.
    fname=@sprintf("q%d.bin",ind[i]);
    q = naive_petsc_read(joinpath(data_dir,fname))/(sqrt(2*pi)*sqrt(dtn_a))
    Q[:,i]= q;
    # Multiply it out (a bit complicated since since we want to save memory)
    # NNZ part of the q-vector
    qnz = q[(start_dtn_dofs+1):end]
    Qnz=sparse(qnz*qnz')
    (I,J,V)=findnz(Qnz)
    P[i]=sparse(I.+start_dtn_dofs,J.+start_dtn_dofs,V);
end



struct DtN_NEP <: NEP
    A::SparseMatrixCSC{ComplexF64}
    M::SparseMatrixCSC{ComplexF64}
    Q::Matrix{ComplexF64};
    P::Vector{SparseMatrixCSC{ComplexF64}} # Full matrix coefficents
    ind2::Vector;
    n::Int;
end



size(nep::DtN_NEP) = (nep.n,nep.n)
size(nep::DtN_NEP,dim::Int) = nep.n


function besselh_quotient(nu,S)
# Compute: besselh'(nu,s)/besselh(nu,s)
    Fder=0.5*(besselh(nu-1,S)-besselh(nu+1,S));
    return besselh(nu,S)\Fder;
end

nep=DtN_NEP(A,M,Q,P,ind2,n);


der=0; λ=1+1im;
function compute_Mder(nep::DtN_NEP,λ::Number,der::Int=0)
    if (der>0)
        error("Not implemented");
    end
    A=copy(nep.A);
    A -= nep.M*λ^2;
    a=1.0;
    for i=1:length(nep.ind2)
        # add:  - s B'_m(s)/B_m(s)
        m=nep.ind2[i];
        f= S ->  -S*besselh_quotient(m,a*S);
        A += f(λ)*nep.P[i]
    end
    return A;
end



# Reference eigvals computed with MATLAB code
λv=[4.370036082655539 - 1.526561686254116im
    2.165520379441414 - 0.537312290190438im
    3.958068685793646 - 0.528959556453958im
    4.397326877370154 - 0.671081974763288im
    1.098616611079115 - 1.005745095689971im]

λ=λv[end];

MM=compute_Mder(nep,λ);
using LinearAlgebra;
z=MM\z; z=z/norm(z);
z=MM\z; z=z/norm(z);
should_be_zero=norm(MM*z)
@show should_be_zero
