export iar_chebyshev
using IterativeSolvers

#TODO: shift-and-scale for PEP
#TODO: test compute_y0 as input

# Types specifying which way to compute y0 in chebyshev iar
abstract type ComputeY0Cheb end;
abstract type ComputeY0ChebDEP <: ComputeY0Cheb end;
abstract type ComputeY0ChebPEP <: ComputeY0Cheb end;
abstract type ComputeY0ChebAuto <: ComputeY0Cheb end;

# Data collected in a precomputation phase
abstract type AbstractPrecomputeData end
type PrecomputeData <: AbstractPrecomputeData
    Tc; L; Ttau; P; P_inv; α; γ; σ
end

"""
    iar_chebyshev(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure=default_errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=1,][orthmethod=DGKS][a=-1,][b=1])

### Infinite Arnoldi method

Runs the infinite Arnoldi method (Chebyshev version) which tries to find eigenvalues close to the shift σ. The shifted-and-scaled Chebyshev polynomials in the interval [a,b] are used as Krylov space.


# Example
```julia-repl
julia> using NonlinearEigenproblems: NEPSolver, NEPCore, Gallery
julia> nep=nep_gallery("dep0");
julia> λ,v=iar_chebyshev(nep);
julia> minimum(svdvals(compute_Mder(nep,λ[1]))) % Is it an eigenvalue?

```

# References
* Algorithm 2 in Jarlebring, Michiels Meerbergen, A linear eigenvalue algorithm for the nonlinear eigenvalue problem, Numer. Math, 2012
"""
iar_chebyshev(nep::NEP;params...)=iar_chebyshev(Complex128,nep;params...)
function iar_chebyshev{T,T_orth<:IterativeSolvers.OrthogonalizationMethod,
                       T_y0<:ComputeY0Cheb}(
    ::Type{T},
    nep::NEP;
    orthmethod::Type{T_orth}=DGKS,
    maxit=30,
    linsolvercreator::Function=default_linsolvercreator,
    tol=eps(real(T))*10000,
    Neig=6,
    errmeasure::Function=default_errmeasure(nep::NEP),
    σ=zero(T),
    γ=one(T),
    v=randn(real(T),size(nep,1)),
    displaylevel=0,
    check_error_every=1,
    compute_y0_method::Type{T_y0}=ComputeY0ChebAuto,
    a=-1.0,
    b=1.0
    )
    #TODO: this function does not work for shift and scaled problems.
    # Possible fix (only for DEP): the shift-and-scaling influences
    # the matrix coefficients and the delay, then we use the usual formulas

    if (compute_y0_method == ComputeY0ChebAuto)
        if (isa(nep,DEP))
             compute_y0_method=ComputeY0ChebDEP;
        elseif (isa(nep,PEP))
            compute_y0_method=ComputeY0ChebPEP;
        else
            compute_y0_method=ComputeY0Cheb;
        end
    end

    #TODO: working on this: explicit shift and scaling
    if ( σ!=zero(T) || γ!=one(T) ) && isa(nep,DEP) || isa(nep,PEP)
            warn("The problem will be explicitly shifted and scaled. The shift and scaling feature is not supported in the general version of iar_chebyshev.")
            # TODO: use the original errmeasure and not compute_resnorm. I don't know why doesn't work
            errmeasure=function (μ,v) return compute_resnorm(nep,σ+γ*μ,v) end
            nep=shift_and_scale(nep,shift=σ,scale=γ);
            σ_orig=σ; γ_orig=γ
            σ=zero(T); γ=one(T)
    end

    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis

    n = size(nep,1);
    m = maxit;

    # initialization
    V = zeros(T,n*(m+1),m+1);
    H = zeros(T,m+1,m);
    y = zeros(T,n,m+1);
    α=γ.^(0:m); α[1]=zero(T);
    local M0inv::LinSolver = linsolvercreator(nep,σ);
    err = ones(m,m);
    λ=zeros(T,m+1); Q=zeros(T,n,m+1);

    vv=view(V,1:1:n,1); # next vector V[:,k+1]
    v=ones(n,1);        # debug
    vv[:]=v; vv[:]=vv[:]/norm(vv);
    k=1; conv_eig=0;

    # hardcoded matrix L
    L=diagm(vcat(2, 1./(2:m)),0)+diagm(-vcat(1./(1:(m-2))),-2); L=L*(b-a)/4;


    # precomputation for exploiting the structure DEP, PEP, GENERAL
    precomp=precompute_data(T,nep,compute_y0_method,a,b,maxit,γ,σ)

    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && (rem(k,check_error_every)==0) || (k==m)
            println("Iteration:",k, " conveig:",conv_eig)
        end
        VV=view(V,1:1:n*(k+1),1:k); # extact subarrays, memory-CPU efficient
        vv=view(V,1:1:n*(k+1),k+1); # next vector V[:,k+1]

        # compute y (only steps different then standard IAR)
        y=zeros(T,n,k+1);
        if (compute_y0_method != ComputeY0Cheb)
            y[:,2:k+1]  = reshape(VV[1:1:n*k,k],n,k)*L[1:k,1:k];
        end
        y[:,1] = compute_y0_cheb(T,nep,compute_y0_method,reshape(VV[1:1:n*k,k],n,k),
                                 y,M0inv,precomp);


        vv[:]=reshape(y[:,1:k+1],(k+1)*n,1);
        # orthogonalization
        H[k+1,k] = orthogonalize_and_normalize!(VV, vv, view(H,1:k,k), orthmethod)

        # compute Ritz pairs (every check_error_every iterations)
        if (rem(k,check_error_every)==0)||(k==m)
            D,Z=eig(H[1:k,1:k]);
            VV=view(V,1:1:n,1:k);
            Q=VV*Z; λ=σ+γ./D;
            conv_eig=0;
            for s=1:k
                err[k,s]=errmeasure(λ[s],Q[:,s]);
                if err[k,s]<tol; conv_eig=conv_eig+1; end
            end
            idx=sortperm(err[k,1:k]); # sort the error
            err[1:k,k]=err[idx,k];
            # extract the converged Ritzpairs
            if (k==m)||(conv_eig>=Neig)
                λ=λ[idx[1:min(length(λ),Neig)]]
                Q=Q[:,idx[1:length(λ)]]
            end
        end
        k=k+1;
    end

    # NoConvergenceException
    if conv_eig<Neig
        err=err[end,1:Neig];
        idx=sortperm(err); # sort the error
        λ=λ[idx];  Q=Q[:,idx]; err=err[idx];
        msg="Number of iterations exceeded. maxit=$(maxit)."
        if conv_eig<3
            msg=string(msg, " Check that σ is not an eigenvalue.")
        end
        throw(NoConvergenceException(λ,Q,err,msg))
    end

    # eventually shift and rescale the eigenvalues if the problem was shifted and rescaled
    try
        λ=σ_orig+γ_orig*λ
    end

    k=k-1
    return λ,Q,err[1:k,:],V[:,1:k]
end

function PrecomputeData()
    return PrecomputeData(0,0,0,0,0,0,0,0);
end

# Precompute data (depending on NEP-type and y0 computation method)
function precompute_data(T,nep::NEPTypes.DEP,::Type{ComputeY0ChebDEP},a,b,m,γ,σ)
    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis
    precomp=PrecomputeData();

    # creating the vector containing T_i(c)
    precomp.Tc=cos.((0:m)'.*acos(cc));  # vector containing T_i(c)

    # creating the matrix containing T_i(-kk*nep.tauv[j]+cc)
    precomp.Ttau=zeros(T,length(nep.tauv),m+2);    II=(0:m+1)';
    for j=1:length(nep.tauv)
        tauv_SS=-kk*nep.tauv[j]+cc;
        if abs(tauv_SS)<1
            precomp.Ttau[j,:]=cos.(II.*acos(tauv_SS));
        elseif tauv_SS>1
            precomp.Ttau[j,:]=cosh.(II.*acosh(tauv_SS));
        else
            precomp.Ttau[j,:]=((-1).^II).*cosh.(II.*acosh(-tauv_SS));
        end
    end
    return precomp;
end
function precompute_data(T,nep::NEPTypes.PEP,::Type{ComputeY0ChebPEP},a,b,m,γ,σ)
    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis
    precomp=PrecomputeData();
    precomp.Tc=cos.((0:m)'.*acos(cc));  # vector containing T_i(c)
    L=diagm(vcat(2, 1./(2:m)),0)+diagm(-vcat(1./(1:(m-2))),-2); L=L*(b-a)/4;
    precomp.L=L

    return precomp;
end
function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb},a,b,m,γ,σ)
    warn("The nep does not belong to the class of DEP or PEP and the function compute_y0 is not provided. Check if the nep belongs to such classes and define it accordingly or provide the function compute_y0. If none of these options are possible, the method will be based on the convertsion between Chebyshev and monomial base and may be numerically unstable if many iterations are performed.")
    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis
    precomp=PrecomputeData();
    L=diagm(vcat(2, 1./(2:m)),0)+diagm(-vcat(1./(1:(m-2))),-2); L=L*(b-a)/4;
    precomp.L=L
    precomp.P=mapslices(x->cheb2mon(T,kk,cc,x),eye(T,m+1,m+1),1)        # P maps chebyshev to monomials as matrix vector action
    precomp.P_inv=mapslices(x->mon2cheb(T,kk,cc,x),eye(T,m+1,m+1),1)    # P_inv maps monomials to chebyshev as matrix vector action
    precomp.σ=σ;
    precomp.γ=γ;
    precomp.α=γ.^(0:m);
    return precomp;
end

function compute_y0_cheb(T,nep::NEPTypes.DEP,::Type{ComputeY0ChebDEP},x,y,M0inv,precomp::AbstractPrecomputeData)
# compute_y0_dep computes y0 for the DEP
# The formula is explicitly given by
# y_0= \sum_{i=1}^N T_{i-1}(γ) x_i - \sum_{j=1}^m A_j \left( \sum_{i=1}^{N+1} T_{i-1}(-ρ \tau_j+γ) y_i\right )
# where T_i is the i-th Chebyshev polynomial of the first kind

    Tc=precomp.Tc;
    Ttau=precomp.Ttau;

    N=size(x,2);   n=size(x,1);
    y0=sum(broadcast(*,x,view(Tc,1:1,1:N)),2); # \sum_{i=1}^N T_{i-1}(γ) x_i
    for j=1:length(nep.tauv) # - \sum_{j=1}^m A_j \left( \sum_{i=1}^{N+1} T_{i-1}(-ρ \tau_j+γ) y_i\right )
        y0-=nep.A[j]*sum(broadcast(*,y,view(Ttau,j:j,1:N+1)),2);
    end
    y0=lin_solve(M0inv,y0)
    return y0
end

function compute_y0_cheb(T,nep::NEPTypes.PEP,::Type{ComputeY0ChebPEP},x,y,M0inv,precomp::AbstractPrecomputeData)
# compute_y0_pep computes y0 for the PEP
# The formula is explicitly given by
# y_0= \sum_{j=0}^{d-1} A_{j+1} x D^j T(c) - y T(c)
# where T(c) is the vector containing T_i(c) as coefficients, where T_i is the i-th Chebyshev polynomial of the first kind

    Tc=precomp.Tc;
    L=precomp.L;
    N=size(x,2);   n=size(x,1);    T=eltype(y);

    # compute the derivation matrix
    LL=inv(L[1:N,1:N])
    D=vcat(zeros(1,N),LL[1:N-1,:]);

    d=length(nep.A)-1;  # degree of the PEP

    # sum for every coefficiet
    v=view(Tc,1:N);
    y0=zeros(T,n,1);
    for j=0:d-1
        y0+=nep.A[j+2]*(x*v);
        v=D*v;
    end
    y0=-lin_solve(M0inv,y0)
    y0-=y*(view(Tc,1:N+1));
    return y0
end

function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb},x,y,M0inv,precomp::AbstractPrecomputeData)
    k=size(x,2);
    α=precomp.α; σ=precomp.σ;
    y[:,2:k+1] = x*precomp.P[1:k,1:k]';
    broadcast!(/,view(y,:,2:k+1),view(y,:,2:k+1),(1:k)')
    y[:,1] = compute_Mlincomb(nep,σ,y[:,1:k+1],a=α[1:k+1]);
    y[:,1] = -lin_solve(M0inv,y[:,1]);
    y[:,:]=y*precomp.P_inv[1:k+1,1:k+1]';
    return y[:,1];
end

function mon2cheb(T,ρ,γ,a)
#cheb2mon: Monomial basis to shifted-and-scaled Chebyshev basis conversion.
#    c = cheb2mon(rho,gamma,a) converts a polynomial written in
#    monomial basis
#    p(x)=\sum_{j=0}^n a(j+1) x^j
#    to a polynomial written in shifted-and-scaled Chebyshev basis
#    p(x)=\sum_{j=0}^n c(j+1) T_j(rho*x+gamma)
#    where T_j(x) is the j-th Chebyshev polynomial of the first kind.
#
#    Example:
#    Suppose we have a polynomial in the monomial basis:
#    2-4*x-x^2
#    This polynomial is expressed in shifted-and-scaled Chebyshev basis as:
#    c(1)*T_0(2*x+3)+c(2)*T_1(2*x+3)+a(3)*T_2(2*x+3)
#    where
#    c=cheb2mon(2,3,[2 -4 -1])
#
#    This function is based on
#    MATH77 and mathc90, Release 6.0, Libraries of Mathematical Subprograms in Fortran 77 and C
#    Fred T. Krogh, Charles L. Lawson, W. Van Snyder
#    in chapter 11 section 3
#    available in the following link
#    http://www.netlib.org/math/docpdf/ch11-03.pdf

    n=length(a)-1;
    α=1/(2*ρ);    β=-γ/ρ;
    b=zeros(T,n+3,1);
    bb=zeros(T,n+3,1);

    for j=n:-1:0
        bb[1]=α*b[2]+β*b[1]+a[j+1];
        bb[2]=β*b[2]+α*b[3]+2*α*b[1];
        for k=3:n-j-1
            bb[k]=α*b[k-1]+β*b[k]+α*b[k+1];
        end
        if n-j>2
            bb[n-j]=α*b[n-j-1]+β*b[n-j];
        end
        if n-j+1>2
            bb[n-j+1]=α*b[n-j];
        end
        b=bb;
        bb=0*bb;
    end
    c=b[1:n+1];
    c=c[:,1];
end

function cheb2mon(T,ρ,γ,c)
#cheb2mon: shifted-and-scaled Chebyshev basis conversion to Monomial basis.
#    a = cheb2mon(T,rho,gamma,c) converts a polynomial written in
#    shifted-and-scaled Chebyshev basis
#    p(x)=\sum_{j=0}^n c(j+1) T_j(rho*x+gamma)
#    to a polynomial written in monomial basis
#    p(x)=\sum_{j=0}^n a(j+1) x^j
#    where T_j(x) is the j-th Chebyshev polynomial of the first kind.
#    T is the type, e.g., Float64
#
#    Example:
#    Suppose we have a polynomial in the shifted-and-scaled Chebyshev basis:
#    2*T_0(2x+3)-4*T_1(2x+3)-T_2(2x+3)
#    This polynomial is expressed in monomial basis as
#    a(1)+a(2)*x+a(3)*x^2
#    where
#    a=cheb2mon(Float64,2,3,[2 -4 -1])
#
#    This function is based on
#    MATH77 and mathc90, Release 6.0, Libraries of Mathematical Subprograms in Fortran 77 and C
#    Fred T. Krogh, Charles L. Lawson, W. Van Snyder
#    in chapter 11 section 3
#    available in the following link
#    http://www.netlib.org/math/docpdf/ch11-03.pdf
    n=length(c)-1;
    α=1/(2*ρ);    β=-γ/ρ;
    a=zeros(T,n+3,1);     b=zeros(T,n+3,1);
    bb=zeros(T,n+3,1);    bb[1:n+1]=c;
    for j=1:1:n+1
        for k=n-j+1:-1:2
            b[k]=(bb[k+1]-β*b[k+1]-α*b[k+2])/α;
        end
        b[1]=(bb[2]-β*b[2]-α*b[3])/(2*α);
        a[j]=bb[1]-α*b[2]-β*b[1];

        bb=b; b=0*b;
    end
    a=a[1:n+1,1];
end
