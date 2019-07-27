
struct ChebPEPFull <: AbstractSPMF{AbstractMatrix}
    n::Int # size
    a; b # Chebyshev polys are scaled to the interval
    Fk::Array   # function values in Chebyshev nodes
    spmf::SPMF_NEP; # The cheb-coefficents are stored directly in the SPMF
end



#type ChebPEPEco <: AbstractSPMF
#    n::Int # size
#    a,b # Chebyshev polys are scaled to the interval
#    Fk::Array   # function values in Chebyshev nodes
#end
#
#ChebPEP = Union{ChebPEPFull}
#ChebPEP = ChebPEPFull;

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

function chebyshev_compute_coefficients(a,b,Fk,k)
    # Return coefficient
    xk=get_chebyshev_nodes(Float64,a,b,k)
    if (size(xk) != size(Fk))
        error("Incompatible sizes");
    end
    Tmat=zeros(k,k);
    for i=1:k
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

# Evaluate a ChebPEP in a point: <=> Interpolate
# Compute coefficients:
# http://inis.jinr.ru/sl/M_Mathematics/MRef_References/Mason,%20Hanscomb.%20Chebyshev%20polynomials%20(2003)/C0355-Ch06.pdf


f=s-> ones(3)-(1:3)*sin(s);
#f=s-> 1-1/(s+6);
#f=s-> 1+s;
k=3;
a=-1;b=4;

(Fk,xk)=chebyshev_eval(a,b,k,f);

Ck0=chebyshev_compute_coefficients(a,b,Fk,k);
#Ck1=chebyshev_compute_coefficients_naive(a,b,Fk,k);

Ck=Ck1;
xk=get_chebyshev_nodes(Float64,a,b,k)
xx=xk[1];
ff0=mapreduce(i->Ck0[i]*cheb_f(a,b,xx,i-1), +, 1:k)
#ff1=mapreduce(i->Ck1[i]*cheb_f(a,b,xx,i-1), +, 1:k)

@show norm(ff0-f(xx))
#@show ff1-f(xx)




function ChebPEP()

end


#function polyeig(pep::ChebPEP)
#
#end

# Compute Mlincomb?


"""

Interpolates the orgspmf in the interval [a,b] with
k chebyshev nodes and gives a representation
in terms of a Chebyshev basis.
"""
function ChebPEPFull(orgspmf,a,b,k)
    F=s-> compute_Mder(orgspmf,s);
    (Fk,xk)=chebyshev_eval(a,b,k,F);
    @show Fk
    @show xk
    Ck=chebyshev_compute_coefficients(a,b,Fk,k);

    fv=Array{Function}(undef,0);
    for j=1:k
        push!(fv, S-> cheb_f(a,b,S,j-1))
    end
    cheb_spmf=SPMF_NEP(Ck,fv);

    n=size(Fk[1],1);
    return ChebPEPFull(n,a,b,Fk,cheb_spmf);
end

nep=nep_gallery("dep0");
nep2=ChebPEPFull(nep,-2,3,23);
nep3=nep2.spmf;

s=0;

compute_Mder(nep,s)-compute_Mder(nep3,s)

λ1=iar(nep,neigs=2,logger=1)[1]
λ3=iar(nep3,neigs=2,logger=1)[1]


@show norm(λ1-λ3)
