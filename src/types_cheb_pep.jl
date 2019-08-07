
export ChebPEP

# Returns the chebyshev nodes scaled to interval [a,b]
# computed with type T
function get_chebyshev_nodes(::Type{T},a,b,k) where {T<:Real}
    mypi=T(pi);
    return (T(a)+T(b))/2 .+ (T(b)-T(a))*(cos.((2*Vector(1:k).-1)*mypi/(2*k)))/2
end

import Base.acos
function acos(S::LowerTriangular)
    # Specialized acos-function for a lower triangular matrix
    # Some extra allocations, but preserves realness.
    #
    # This is a work-around julia issue:
    #    https://github.com/JuliaLang/julia/issues/32721
    F=acos(Matrix(S));
    if (all(isreal.(acos.(diag(S)))))
        F .= real(F);
    end
end
function cheb_f_cosine_formula(a,b,x,k)
    x1=2*(x-a*one(x))/(b-a)-one(x);
    if (istril(x1) && (x1 isa AbstractMatrix))
        x1=LowerTriangular(x1);
    end
    return cos(k*acos(x1));
end
function cheb_f_poly(a,b,x,k)
    x1=2*(x-a*one(x))/(b-a) -one(x);
    if (k==0)
        return one(x1);
    elseif k==1
        return x1
    elseif k==2
        return 2*x1^2-one(x1)
    elseif k==3
        return 4*x1^3-3*x1
    elseif k==4
        return 8*x1^4-8*x1^2+one(x1)
    elseif k==5
        return 16*x1^5 - 20*x1^3 + 5*x1;
    elseif k==6
        return 32*x1^6-48*x1^4+18*x1^2-one(x1);
    elseif k==7
        return 64*x1^7-112*x1^5+56*x1^3-7*x1;
    elseif k==8
        return 128*x1^8-256*x1^6+160*x1^4-32*x1^2+one(x1);
    else
        error("Not implemented")
    end
end

function cheb_f(a,b,x,k)
    return cheb_f_cosine_formula(a,b,x,k);
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
# Chebyshev Polynomials, 1st Edition, J.C. Mason, David C. Handscomb
# Chapter 8
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
    ChebPEP(orgnep::NEP,k,[a=-1,[b=1]] [,cosine_formula_cutoff=5])

The type `ChebPEP<:AbstractSPMF` represents a polynomial function
where the
function is stored using a Chebyshev basis scaled to the
interval `[a,b]`, i.e.,
```math
M(λ)= B_0T_0(λ)+⋯+B_{k-1}T_{k-1}(λ)
```
where `T_i` are the scaled and shifted Chebyshev polynomials.

The creator `ChebPEP` takes `nep::NEP` as an input
and interpolates this NEP in `k` Chebyshev nodes, resulting
in a polynomial of degree `k-1`, represented by its
coefficients in the Chebyshev basis.
Interpolation in Chebyshev nodes is known to have
attractive approximation properties, as well
as robustness with respect to round-off errors.

The kwarg `cosine_formula_cutoff` decides how the Chebyshev
polynomials should be computed. For larger degrees, one
wants to use the cosine formula, whereas for low degrees
the explicit monomial expression is more efficient.
The explicit monomial expression will be used for degrees
lower than `cosine_formula_cutoff`.


# Example:

```julia
julia> nep=nep_gallery("dep0");
julia> chebpep=ChebPEP(nep,9);
julia> using LinearAlgebra;
julia> norm(compute_Mder(nep,0.3)-compute_Mder(chebpep,0.3))
1.6920305863798614e-8
julia> chebpep=ChebPEP(nep,19); # Better interpolation
julia> norm(compute_Mder(nep,0.3)-compute_Mder(chebpep,0.3))
2.2782119183158333e-15
```

See also: `polyeig`
"""
function ChebPEP(orgspmf,k,a=-1,b=1;cosine_formula_cutoff=5)
    F=s-> compute_Mder(orgspmf,s);
    (Fk,xk)=chebyshev_eval(a,b,k,F);
    Ck=chebyshev_compute_coefficients(a,b,Fk,xk);

    fv=Array{Function}(undef,0);
    for j=1:k
        if (j>cosine_formula_cutoff)
            push!(fv, S-> cheb_f_cosine_formula(a,b,S,j-1))
        else
            push!(fv, S-> cheb_f_poly(a,b,S,j-1))
        end
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
