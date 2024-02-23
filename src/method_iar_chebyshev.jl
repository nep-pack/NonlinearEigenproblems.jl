export iar_chebyshev
export compute_y0_cheb

using IterativeSolvers
using LinearAlgebra
using Random

# Types specifying which way to compute y0 in chebyshev iar
abstract type ComputeY0Cheb end;
abstract type ComputeY0ChebDEP <: ComputeY0Cheb end;
abstract type ComputeY0ChebPEP <: ComputeY0Cheb end;
abstract type ComputeY0ChebSPMF_NEP <: ComputeY0Cheb end;
abstract type ComputeY0ChebAuto <: ComputeY0Cheb end;


# Data collected in a precomputation phase.
# These are made mutable (could be made immutable by appropriate modification in precompute_data)
abstract type AbstractPrecomputeData end
mutable struct PrecomputeDataDEP <: AbstractPrecomputeData
    Tc; Ttau;
end
mutable struct PrecomputeDataPEP <: AbstractPrecomputeData
    Tc; D;
end
mutable struct PrecomputeDataSPMF <: AbstractPrecomputeData
    Tc; DDf;
end
mutable struct PrecomputeDataNEP <: AbstractPrecomputeData
    P; P_inv; α; γ; σ;
end


"""
    iar_chebyshev(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=DefaultLinSolverCreator(),][tolerance=eps()*10000,][neigs=6,][errmeasure,][v=rand(size(nep,1),1),][logger=0,][check_error_every=1,][orthmethod=DGKS][a=-1,][b=1,][compute_y0_method=ComputeY0ChebAuto])

Run the infinite Arnoldi method (Chebyshev version) on the nonlinear eigenvalue problem stored in `nep`.

The target `σ` is the center around which eiganvalues are computed. A Ritz pair `λ` and `v` is flagged a as converged (to an eigenpair) if `errmeasure` is less than `tol`. The vector
`v` is the starting vector for constructing the Krylov space. The orthogonalization method, used in contructing the orthogonal basis of the Krylov space, is specified by `orthmethod`, see the package `IterativeSolvers.jl`. The iteration
is continued until `neigs` Ritz pairs converge. This function throws a `NoConvergenceException` if the wanted eigenpairs are not computed after `maxit` iterations.
However, if `neigs` is set to `Inf` the iteration is continued until `maxit` iterations without an error being thrown.
The kwarg `compute_y0_method` specifying how the next vector of the Krylov space (in Chebyshev format) can be computed. See [`compute_y0_cheb`](@ref) in the module NEPSolver with the command `?NEPSolver.compute_y0_cheb`.

See [`augnewton`](@ref) for other parameters.


# Example
```julia-repl
julia> using NonlinearEigenproblems, LinearAlgebra
julia> nep=nep_gallery("dep0",100);
julia> v0=ones(size(nep,1));
julia> λ,v=iar_chebyshev(nep;v=v0,tol=1e-5,neigs=3);
julia> norm(compute_Mlincomb!(nep,λ[1],v[:,1])) # Is it an eigenvalue?
julia> λ    # print the computed eigenvalues
3-element Array{Complex{Float64},1}:
julia> norm(compute_Mlincomb(nep,λ[1],v[:,1]))
 0.050462487848960284 - 1.4289626573515395e-18im
 -0.07708779190301127 + 7.703053374113074e-18im
   0.1503856540695659 - 1.662582577182149e-17im
```

# References
* Algorithm 2 in Jarlebring, Michiels Meerbergen, A linear eigenvalue algorithm for the nonlinear eigenvalue problem, Numer. Math, 2012
"""
iar_chebyshev(nep::NEP;params...)=iar_chebyshev(ComplexF64,nep;params...)
function iar_chebyshev(
    ::Type{T},
    nep::NEP;
    orthmethod=DGKS(),
    maxit=30,
    linsolvercreator=DefaultLinSolverCreator(),
    tol=eps(real(T))*10000,
    neigs=6,
    errmeasure::ErrmeasureType = DefaultErrmeasure(nep),
    σ=zero(T),
    γ=one(T),
    v=randn(real(T),size(nep,1)),
    logger=0,
    check_error_every=1,
    compute_y0_method::Type{T_y0}=ComputeY0ChebAuto,
    a = isa(nep,DEP) ? -maximum(nep.tauv) : -1.0,
    b = isa(nep,DEP) ? 0.0 : 1.0
    )where{T,T_y0<:ComputeY0Cheb}

    @parse_logger_param!(logger)

    if (compute_y0_method == ComputeY0ChebAuto)
        if (isa(nep,DEP))
            compute_y0_method=ComputeY0ChebDEP;
        elseif (isa(nep,PEP))
            compute_y0_method=ComputeY0ChebPEP;
        elseif  (isa(nep,SPMF_NEP))
            compute_y0_method=ComputeY0ChebSPMF_NEP;
        else
            compute_y0_method=ComputeY0Cheb;
        end
    end

    if ( σ!=zero(T) || γ!=one(T) ) && compute_y0_method<:Union{ComputeY0ChebDEP,ComputeY0ChebPEP}
            @warn "The problem will be explicitly shifted and scaled. The shift and scaling feature is not supported in the general version of iar_chebyshev."
            # TODO: use the original errmeasure and not compute_resnorm. I don't know why doesn't work
            #errmeasure=(μ,v)-> errmeasure(σ+γ*μ,v) #this is what we want but it does not work
            errmeasure=function (μ,v) compute_resnorm(nep,σ+γ*μ,v) end
            nep=shift_and_scale(nep,shift=σ,scale=γ);
            σ_orig=σ; γ_orig=γ
            σ=zero(T); γ=one(T)
    end
    push_info!(logger,"IAR Chebyshev with interval [$a,$b]");

    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis

    n = size(nep,1);
    m = maxit;


    # initialization
    V = zeros(T,n*(m+1),m+1);
    H = zeros(T,m+1,m);
    y = zeros(T,n,m+1);
    α=Vector{T}(γ.^(0:m)); α[1]=zero(T);
    local M0inv::LinSolver=create_linsolver(linsolvercreator,nep,σ)
    err = ones(m,m);
    λ=zeros(T,m+1); Q=zeros(T,n,m+1);

    vv=view(V,1:1:n,1); # next vector V[:,k+1]
    vv[:]=v; vv[:]=vv[:]/norm(vv);
    k=1; conv_eig=0;

    # hardcoded matrix L
    L=diagm(0 => vcat(2, 1 ./ (2:m))) + diagm(-2 => -vcat(1 ./ (1:(m-2))))
    L=L*(b-a)/4

    # precomputation for exploiting the structure DEP, PEP, GENERAL
    precomp=precompute_data(T,nep,compute_y0_method,a,b,maxit,γ,σ)


    while (k <= m) && (conv_eig<neigs)

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
        if ((rem(k,check_error_every)==0)||(k==m))&&(k>2)
            D,Z = eigen(H[1:k,1:k])
            VV=view(V,1:1:n,1:k);
            Q=VV*Z
            λ=σ .+ γ ./ D
            conv_eig=0;
            # compute the errors
            err[k,1:size(λ,1)]=
              map(s-> estimate_error(errmeasure,λ[s],Q[:,s]), 1:size(λ,1))
            # Log them and compute the converged
            push_iteration_info!(logger,2, k,err=err[k,1:size(λ,1)],
                                 continues=true);
            for s=1:size(λ,1)
                if err[k,s]<tol;
                    conv_eig=conv_eig+1;
                    push_info!(logger,"+", continues=true);
                elseif err[k,s]<tol*10
                    push_info!(logger,"=", continues=true);
                else
                    push_info!(logger,"-", continues=true);
                end
            end
            push_info!(logger,"");
            # Sort the errors
            idx=sortperm(err[k,1:k]); # sort the error
            err[k,1:k]=err[k,idx];

            # extract the converged Ritzpairs
            if (k==m)||(conv_eig>=neigs)
                nrof_eigs = Int(min(length(λ),neigs))
                λ=λ[idx[1:nrof_eigs]]
                Q=Q[:,idx[1:nrof_eigs]]
            end
        end
        k=k+1;
    end
    k=k-1
    # NoConvergenceException
    if conv_eig<neigs && neigs != Inf
        err=err[end,1:neigs];
        idx=sortperm(err); # sort the error
        λ=λ[idx];  Q=Q[:,idx]; err=err[idx];
        msg="Number of iterations exceeded. maxit=$(maxit)."
        if conv_eig<3
            msg=string(msg, " Check that σ is not an eigenvalue.")
        end
        throw(NoConvergenceException(λ,Q,err,msg))
    end

    # eventually shift and rescale the eigenvalues if the problem was shifted and rescaled
    # TODO: we shouldn't use the exception system for this, perhaps a Bool flag is a better choice
    try
        λ = σ_orig .+ γ_orig * λ
    catch
        # ignore
    end

    # extract the converged Ritzpairs
    λ=λ[1:min(length(λ),conv_eig)];
    Q=Q[:,1:min(size(Q,2),conv_eig)];
    return λ,Q,err[1:k,:],V[:,1:k],H[1:k,1:k]
end

# Contructors for the precomputed data
function PrecomputeDataInit(::Type{ComputeY0ChebDEP})
    return PrecomputeDataDEP(0,0);
end
function PrecomputeDataInit(::Type{ComputeY0ChebPEP})
    return PrecomputeDataPEP(0,0);
end
function PrecomputeDataInit(::Type{ComputeY0ChebSPMF_NEP})
    return PrecomputeDataSPMF(0,0);
end
function PrecomputeDataInit(::Type{ComputeY0Cheb})
    return PrecomputeDataNEP(0,0,0,0,0);
end

# Precompute data (depending on NEP-type and y0 computation method)
function precompute_data(T,nep::NEPTypes.DEP,::Type{ComputeY0ChebDEP},a,b,m,γ,σ)
    if ( σ!=zero(T) || γ!=one(T) )
        error("This function does not support shift and scale parameters");
    end
    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis
    precomp=PrecomputeDataInit(ComputeY0ChebDEP);
    # creating the vector containing T_i(c)
    precomp.Tc=cos.((0:m)'.*acos(cc));  # vector containing T_i(c)
    # creating the matri  containing T_i(-kk*nep.tauv[j]+cc)
    precomp.Ttau=zeros(T,length(nep.tauv),m+2);    II=(0:m+1)';
    for j=1:length(nep.tauv)
        tauv_SS=-kk*nep.tauv[j]+cc;
        if abs(tauv_SS)<=1
            precomp.Ttau[j,:]=cos.(II.*acos(tauv_SS));
        elseif tauv_SS>=1
            precomp.Ttau[j,:]=cosh.(II.*acosh(tauv_SS));
        else
            precomp.Ttau[j,:]=((-1).^II).*cosh.(II.*acosh(-tauv_SS));
        end
    end
    return precomp;
end
function precompute_data(T,nep::NEPTypes.PEP,::Type{ComputeY0ChebPEP},a,b,m,γ,σ)
    # The matrix D is the derivation map, namely D_N [T_0(s),T_1(s), ..., T_{N-1}(s)]=[T_0'(s),T_1'(s), ..., T_{N-1}'(s)]
    # NOTE: The matrix D_N is the submatrix D_M[1:N+1,1:N+1] for M>=N (larger derivarive matrix)
    if ( σ!=zero(T) || γ!=one(T) )
        error("This function does not support shift and scale parameters");
    end
    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis
    precomp=PrecomputeDataInit(ComputeY0ChebPEP);
    precomp.Tc=cos.((0:m)'.*acos(cc));  # vector containing T_i(c)
    L=diagm(0 => vcat(2, 1 ./ (2:m))) + diagm(-2 => -vcat(1 ./ (1:(m-2))))
    L=L*(b-a)/4
    L=inv(L[1:m,1:m]); precomp.D=vcat(zeros(1,m),L[1:m-1,:]);
    return precomp;
end
function precompute_data(T,nep::NEPTypes.AbstractSPMF,::Type{ComputeY0ChebSPMF_NEP},a,b,m,γ,σ)
    # Precomputation for NEPS in SPMF format
    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis
    precomp=PrecomputeDataInit(ComputeY0ChebSPMF_NEP);
    precomp.Tc=cos.((0:m)'.*acos(cc));  # vector containing T_i(c)
    L=diagm(0 => vcat(2, 1 ./ (2:m))) + diagm(-2 => -vcat(1 ./ (1:(m-2))))
    L=L*(b-a)/4; L=inv(L[1:m,1:m]); # integration matrix
    D=vcat(zeros(1,m),L[1:m-1,:]); # detivation matrix

    fv, Av = get_fv(nep), get_Av(nep)
    DDf = Array{Array{T,2}}(undef, length(fv))
    for i = 1:length(fv)
        DDs = Matrix{T}(σ*I, size(D,1), size(D,1)) + γ*D
        DDf[i] = γ * DD0_mat_fun(T, fv[i], DDs, σ)
    end
    precomp.DDf = DDf
    return precomp
end
function precompute_data(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb},a,b,m,γ,σ)
    # Precomputation for a NEP that is not a DEP, PEP or SPMF
    @warn "The nep does not belong to the class of DEP, PEP or SPMF and the function compute_y0 is not provided. Check if the nep belongs to such classes and define it accordingly or provide the function compute_y0. If none of these options are possible, the method will be based on the convertsion between Chebyshev and monomial base and may be numerically unstable if many iterations are performed."
    cc=(a+b)/(a-b);   kk=2/(b-a); # scale and shift parameters for the Chebyshev basis
    precomp=PrecomputeDataInit(ComputeY0Cheb);
    precomp.P = mapslices(x->cheb2mon(T,kk,cc,x), Matrix{T}(I,m+1,m+1), dims = 1)'        # P maps chebyshev to monomials as matrix vector action
    precomp.P_inv = mapslices(x->mon2cheb(T,kk,cc,x), Matrix{T}(I,m+1,m+1), dims = 1)'    # P_inv maps monomials to chebyshev as matrix vector action
    precomp.σ=σ;
    precomp.γ=γ;
    precomp.α=γ.^(0:m);
    return precomp;
end
"""
    y0 = compute_y0_cheb([eltype],nep::NEPTypes.DEP,::Type{ComputeY0ChebPEP},X,Y,M0inv,precomp::AbstractPrecomputeData)

Computes the vector y0 used in [`iar_chebyshev`](@ref) given by
```math
 y_0 = \\sum_{i=1}^N T_{i-1}(γ) x_i - \\sum_{j=1}^m A_j \\left( \\sum_{i=1}^{N+1} T_{i-1}(-ρ \\tau_j+γ) y_i \\right )
```
where T(c) is the vector containing \$T_i(c)\$ as coefficients, where \$T_i\$ is the i-th Chebyshev polynomial of the first kind.
"""
function compute_y0_cheb(T,nep::NEPTypes.DEP,::Type{ComputeY0ChebDEP},X,Y,M0inv,precomp::AbstractPrecomputeData)
    Tc=precomp.Tc;
    Ttau=precomp.Ttau;
    Av=get_Av(nep);

    n,N=size(X)
    y0=sum(broadcast(*,X,view(Tc,1:1,1:N)), dims = 2) # \sum_{i=1}^N T_{i-1}(γ) x_i
    for j=1:length(nep.tauv) # - \sum_{j=1}^m A_j \left( \sum_{i=1}^{N+1} T_{i-1}(-ρ \tau_j+γ) y_i\right )
        y0-=Av[j+1]*sum(broadcast(*,Y,view(Ttau,j:j,1:N+1)), dims = 2)
    end
    y0=lin_solve(M0inv,y0)
    return y0
end
"""
    y0 = compute_y0_cheb([eltype],nep::NEPTypes.PEP,::Type{ComputeY0ChebPEP},X,Y,M0inv,precomp::AbstractPrecomputeData)

Computes the vector y0 used in [`iar_chebyshev`](@ref) given by
```math
 y_0 = \\sum_{j=0}^{d-1} A_{j+1} x D^j T(c) - y T(c)
```
where T(c) is the vector containing \$T_i(c)\$ as coefficients, where \$T_i\$ is the i-th Chebyshev polynomial of the first kind and \$D\$ is the derivation matrix in Chebyshev basis.
"""
function compute_y0_cheb(T,nep::NEPTypes.PEP,::Type{ComputeY0ChebPEP},X,Y,M0inv,precomp::AbstractPrecomputeData)
    Tc=precomp.Tc;
    n,N=size(X)
    d=length(nep.A)-1;  # degree of the PEP
    # sum for every coefficiet
    v=view(Tc,1:N);
    y0=zeros(T,n,1);
    for j=0:d-1
        y0+=nep.A[j+2]*(X*v);
        v=view(precomp.D,1:N,1:N)*v;
    end
    y0=-lin_solve(M0inv,y0)
    y0-=Y*(view(Tc,1:N+1));
    return y0
end
"""
    y0 = compute_y0_cheb([eltype],nep::NEPTypes.SPMF_NEP,::Type{ComputeY0ChebPEP},X,Y,M0inv,precomp::AbstractPrecomputeData)

Computes the vector y0 used in [`iar_chebyshev`](@ref) given by
```math
 y_0= \\sum_{j=0}^{m} M^{(j)}(\\mu) X b_j \\left( D_N \\right) T_N(c) - Y T_N(c)
```
where T(c) is the vector containing \$T_i(c)\$ as coefficients, where \$T_i\$ is the i-th Chebyshev polynomial of the first kind and \$b_j(\\lambda)=(f_j(0)-f_j(\\lambda))/\\lambda=f[\\lambda,0]\$ are divided differences.
"""
function compute_y0_cheb(T,nep::NEPTypes.AbstractSPMF,::Type{ComputeY0ChebSPMF_NEP},X,Y,M0inv,precomp::AbstractPrecomputeData)
    Tc=precomp.Tc;
    n,N=size(X);
    fv,Av=get_fv(nep),get_Av(nep)
    y0=zero(X)
    for i=1:length(fv)
        y0+=Av[i]*X*view(precomp.DDf[i],1:N,1:N)
    end
    y0=y0*(vec(view(Tc,1:1,1:N)))
    y0=-lin_solve(M0inv,y0)
    y0-=Y*(view(Tc,1:N+1));
    return y0
end
"""
    y0 = compute_y0_cheb([eltype],nep::NEPTypes.NEP,::Type{ComputeY0ChebNEP},X,Y,M0inv,precomp::AbstractPrecomputeData)

Computes the vector y0 used in [`iar_chebyshev`](@ref) defined as
```math
 y_0 =\\left( \\sum_{i=0}^{N-1} B \\left( \\frac{d}{d \\theta} \\right) \\hat T_i(\\theta) x_i \\right)(0) - \\sum_{i=0}^{N} T_i(c) y_i
```
where \$T_i\$ is the i-th Chebyshev polynomial of the first kind, \$ \\hat T_i\$ is the i-th Chebyshev polynomial of the first kind for the interval [a,b]. For a generic `nep`, this quantity is computed by converting polynomials in monomial basis. This procedure may be numerical unstable if many iterations are required. If for the specific `nep` a closed formula is available, we suggest to overload this function.
"""
function compute_y0_cheb(T,nep::NEPTypes.NEP,::Type{ComputeY0Cheb},X,Y,M0inv,precomp::AbstractPrecomputeData)
    # compute_y0_cheb computes y0 for the NEP
    # This function convert the coefficients x=vec(X) from the Chebyshev basis to the monomial basis, then compute the vector y=vec(Y) as in iar (Taylor version) and then convert y=vec(Y) in Chebyshev basis
    k=size(X,2);
    α=precomp.α; σ=precomp.σ;
    Y[:,2:k+1] = X*view(precomp.P,1:k,1:k);
    broadcast!(/,view(Y,:,2:k+1),view(Y,:,2:k+1),(1:k)')
    Y[:,1] = compute_Mlincomb(nep,σ,view(Y,:,1:k+1),α[1:k+1]);
    Y[:,1] = -lin_solve(M0inv,Y[:,1]);
    Y[:,:]=Y*view(precomp.P_inv,1:k+1,1:k+1);
    return Y[:,1];
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
function DD0_mat_fun(T,f,S,σ)
    # evaluete the divided differences matrix function
    # f[S,0] by using the equality
    #
    # f(S I) = (f(S) 		f[S,σI]  )
    #  (0 0)   (0			f(0)	)
    #
    # or, written in a more compact way, it is
    #
    # f([S I; 0 σI])=[f(S) f[S,σI]; 0 f(σ)I]
    #
    # Notice that f[S,0] is defined also for S singular.
    # If S is nonsingular it holds f[S,0]=S^(-1)-(f(S)-f(0))
    # Example:
    # n=10; S=rand(n,n); T=ComplexF64; f=x->exp(x)+x^2
    # Y1=DD0_mat_fun(T,f,S); Y2=inv(S)*(f(S)-f(zeros(S)));
    # opnorm(Y1-Y2)

    n=size(S,1);
    A=zeros(T,2*n,2*n);
    A[1:n, 1:n] = S
    A[1:n, n+1:2*n] = Matrix{T}(I, n, n)
    A[n.+(1:n), n.+(1:n)] = Matrix{T}(σ*I, n, n)
    return f(A)[1:n,n+1:end];
end
