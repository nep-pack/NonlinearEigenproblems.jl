
export ChebPEP


# Returns the chebyshev nodes scaled to interval [a,b]
# computed with type T
function get_chebyshev_nodes(::Type{T},a,b,k) where {T<:Real}
    mypi=T(pi);
    return (T(a)+T(b))/2 .+ (T(b)-T(a))*(cos.((2*Vector(1:k).-1)*mypi/(2*k)))/2
end
function cheb_f_cosine_formula(a,b,x,k)
    x1=2*(x-a*one(x))/(b-a)-one(x);
    return cos(k*acos(x1));
end
function cheb_f_poly(a,b,x,k)
    x1=2*(x-a*one(x))/(b-a) -one(x);
    if (k==0)
        return 1;
    elseif k==1
        return x1
    elseif k==2
        return 2*x1^2-one(x1)
    elseif k==3
        return 4*x1^3-3*x1
    elseif k==4
        return 8*x1^4-8*x^2+one(x1)
    else
        error("Not implemented")
    end
end
function cheb_f(a,b,x,k)
    return cheb_f_cosine_formula(a,b,x,k);
#    return cheb_f_poly(a,b,x,k);
end

# Evaluate F in the chebyshev nodes
function chebyshev_eval(a,b,k,F::Function)
    x=get_chebyshev_nodes(Float64,a,b,k)
    Fk=F.(x);
    return (Fk,x)
end


function chebyshev_compute_coefficients_naive(a,b,Fk,k)
    xk=get_chebyshev_nodes(Float64,a,b,k)
    @show xk
    if (size(xk) != size(Fk))
        error("Incompatible sizes");
    end


    Tinv=zeros(k,k);
    for i=1:k
        for j=1:k
            Tinv[i,j]=cheb_f(a,b,xk[i],j-1);
        end
    end
    @show Tinv

    Ck=Tinv\Fk

    return Ck

end

# Compute the chebyshev coefficients of the coefficients
# stored in Fk. xk should be the chebyshev points
function chebyshev_compute_coefficients(a,b,Fk,xk)
    # Return coefficient
    k=size(Fk,1);
    Tmat=zeros(k,k);
    for i=1:k
        # for each chebyshev polynomial, compute the k coefficients
        Tmat[i,:]= map(x->cheb_f(a,b,x,i-1)*2/k,xk)
    end
    Tmat[1,:] *= 0.5; # Fix the prime in the sum

    # Compute the coefficients from the Tmat via
    # a "matrix vector multiplication"

    # When Fk are scalars
    #Ck=Tmat*Fk;
    # Generalization when Fk are vectors or matrices:
    Ck=map(i-> sum(Fk .* Tmat[i,:]), 1:k)


    return Ck
end







struct ChebPEP{T<:AbstractMatrix,Ftype} <: AbstractSPMF{T}
    n::Int # size
    # Chebyshev polys are scaled to the interval [a,b]
    a::Ftype
    b::Ftype
    k::Int; # Number of Chebyshev polys
    spmf::SPMF_NEP{T,Ftype}; # The cheb-coefficents are stored directly in the SPMF
end



# Delegate compute functions
#import NonlinearEigenproblems.NEPCore.compute_Mder
#import NonlinearEigenproblems.NEPCore.compute_Mlincomb
#import NonlinearEigenproblems.NEPCore.compute_MM
#import NonlinearEigenproblems.NEPTypes.get_Av
#import NonlinearEigenproblems.NEPTypes.get_fv
#import NonlinearEigenproblems.NEPSolver.polyeig
#import Base.size
#using LinearAlgebra
#

function size(nep::ChebPEP)
    return (nep.n,nep.n)
end
function size(nep::ChebPEP,dim)
    return nep.n
end

compute_Mlincomb(nep::ChebPEP,λ::Number,V::AbstractVecOrMat,a::Vector)= compute_Mlincomb(nep.spmf,λ,V,a);
compute_Mlincomb(nep::ChebPEP,λ::Number,V::AbstractVecOrMat)= compute_Mlincomb(nep.spmf,λ,V);
compute_Mder(nep::ChebPEP,λ::Number)=compute_Mder(nep.spmf,λ,0)
compute_Mder(nep::ChebPEP,λ::Number,der::Integer)=compute_Mder(nep.spmf,λ,der)
compute_MM(nep::ChebPEP,par...)=compute_MM(nep.spmf,par...)
get_Av(nep::ChebPEP)=get_Av(nep.spmf)
get_fv(nep::ChebPEP)=get_fv(nep.spmf)



"""

Interpolates the orgspmf in the interval [a,b] with
k chebyshev nodes and gives a representation
in terms of a Chebyshev basis.
"""
function ChebPEP(orgspmf,k,a=-1,b=1)
    F=s-> compute_Mder(orgspmf,s);
    (Fk,xk)=chebyshev_eval(a,b,k,F);
    Ck=chebyshev_compute_coefficients(a,b,Fk,xk);

    fv=Array{Function}(undef,0);
    for j=1:k
        push!(fv, S-> cheb_f(a,b,S,j-1))
    end

    # Determine the types (based on Ck[1])
    MatType=typeof(Ck[1]);
    T=eltype(Ck[1]);

    # Create an SPMF with the correct type
    cheb_spmf=SPMF_NEP(Ck,fv,check_consistency=false,Ftype=T);

    n=size(Fk[1],1);
    # Instantiate the type
    return ChebPEP{MatType,T}(n,a,b,k,cheb_spmf);
end


# TODO: Implement interpolation similar to Effenberger and Kressner. "Chebyshev interpolation for nonlinear eigenvalue problems." BIT Numerical Mathematics 52.4 (2012): 933-951.
function polyeig(pep::ChebPEP{T,Ftype}) where {T,Ftype}
    k=pep.k
    Fk=get_Av(pep);
    @show k
    n=size(pep,1);
    L0=zeros(Ftype,n*(k-1),n*(k-1));
    L1=zeros(Ftype,n*(k-1),n*(k-1));
    II=Matrix{Ftype}(I,n,n);
    @show L0
    for j=1:(k-2)
        L0[((j-1)*n) .+ (1:n), j*n .+ (1:n)]=II
        L0[j*n .+ (1:n), ((j-1)*n) .+ (1:n)]=II
    end

    for j=1:k-1
        L0[((k-2)*n) .+ (1:n), (n*(j-1)).+ (1:n)]=-Fk[j];
    end

    @show L0
    L0[((k-2)*n) .+ (1:n), (n*(k-3)) .+ (1:n)] += Fk[k];

    @show L0
    for j=1:k-2
        factor=2;
        if (j==1)
            factor=1
        end
        L1[((j-1)*n) .+ (1:n), ((j-1)*n) .+ (1:n)]=factor*II
    end
    L1[((k-2)*n) .+ (1:n),
       ((k-2)*n) .+ (1:n)]=2*Fk[k]

    LL=eigen(L0,L1);

    @show L1
    return (LL.values,
            hcat(map(i -> normalize(LL.vectors[1:n,i]), 1:(n*(k-1)))...))

end
