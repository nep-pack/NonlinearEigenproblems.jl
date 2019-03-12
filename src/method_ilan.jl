export ilan

using IterativeSolvers
using LinearAlgebra
using Random
using Statistics

# Types specifying which way to B_prod in chebyshev iar
abstract type Compute_Bmul_method end;
abstract type Compute_Bmul_method_DEP <: Compute_Bmul_method end;
abstract type Compute_Bmul_method_SPMF_NEP <: Compute_Bmul_method end;
abstract type Compute_Bmul_method_DerSPMF <: Compute_Bmul_method end;
abstract type Compute_Bmul_method_Auto <: Compute_Bmul_method end;


# Data collected in a precomputation phase.
abstract type IlanAbstractPrecomputeData end
mutable struct IlanPrecomputeDataDEP <: IlanAbstractPrecomputeData
     G; vv; ZZ; QQ; QQ2
end
mutable struct IlanPrecomputeDataSPMF_NEP <: IlanAbstractPrecomputeData
     G; FDH; ZZ; QQ
end
mutable struct IlanPrecomputeDataDerSPMF <: IlanAbstractPrecomputeData
     G; FDH; ZZ; QQ
end

"""
    ilan(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=30,][orthmethod=DGKS])

Run the infinite Lanczos method on the symmetric nonlinear eigenvalue problem stored in `nep`.

The target `σ` is the center around which eiganvalues are computed. The kwarg `errmeasure` is a function handle which can be used to specify how the error is measured to be used in termination (default is absolute residual norm). A Ritz pair `λ` and `v` is flagged a as converged (to an eigenpair) if `errmeasure` is less than `tol`. The vector
`v` is the starting vector for constructing the Krylov space. The orthogonalization method, used in contructing the orthogonal basis of the Krylov space, is specified by `orthmethod`, see the package `IterativeSolvers.jl`. The iteration
is continued until `Neig` Ritz pairs converge. This function throws a `NoConvergenceException` if the wanted eigenpairs are not computed after `maxit` iterations. The `linsolvercreator` is a function which specifies how the linear system is created and solved.


# Example
```julia-repl
julia> using NonlinearEigenproblems, LinearAlgebra
julia> nep=nep_gallery("dep_symm_double",10);
julia> v0=ones(size(nep,1));
julia> λ,v=ilan(nep;v=v0,tol=1e-5,Neig=3);
julia> norm(compute_Mlincomb!(nep,λ[1],v[:,1])) # Is it an eigenvalue?
julia> λ    # print the computed eigenvalues
3-element Array{Complex{Float64},1}:
 0.04103537900075572 - 1.6342212662372832e-19im
 0.04103537900077957 - 2.5916996904875994e-19im
 0.04114919035623714 - 7.9738202040662040e-20im
```

# References
* Algorithm 2 in Mele, The infinite Lanczos method for symmetric nonlinear eigenvalue problems, https://arxiv.org/abs/1812.07557, 2018
"""
ilan(nep::NEP;params...)=ilan(ComplexF64,nep;params...)
function ilan(
    ::Type{T},
    nep::NEP;
    orthmethod::Type{T_orth}=DGKS,
    maxit=30,
    linsolvercreator::Function=default_linsolvercreator,
    tol=eps(real(T))*10000,
    Neig=6,
    errmeasure::ErrmeasureType = DefaultErrmeasure,
    σ=zero(T),
    γ=one(T),
    v=randn(real(T),size(nep,1)),
    displaylevel=0,
    check_error_every=30,
    inner_solver_method=DefaultInnerSolver,
    Compute_Bmul_method::Type{T_y0}=Compute_Bmul_method_Auto,
    )where{T<:Number,T_orth<:IterativeSolvers.OrthogonalizationMethod,T_y0<:Compute_Bmul_method}

    # Ensure types σ and v are of type T
    σ=T(σ)
    v=Array{T,1}(v)

    if (Compute_Bmul_method == Compute_Bmul_method_Auto)
        if (isa(nep,DEP))
            Compute_Bmul_method=Compute_Bmul_method_DEP
        elseif (isa(nep,DerSPMF))
            Compute_Bmul_method=Compute_Bmul_method_DerSPMF
        elseif (isa(nep,SPMF_NEP))
            Compute_Bmul_method=Compute_Bmul_method_SPMF_NEP
        else
            println("Error in Compute_Bmul_method. You provided and invalid `nep` for which `Compute_Bmul_method` is not automatically supported. Provide the `nep` in a compatible format: `DEP`, `SPMF_NEP` or `DerSPMF`. Otherwise you can overload the function `Bmult` for your specif `nep`.")
        end
    end


    n = size(nep,1);
    m = maxit;

    # initialization
    V=zeros(T,n,m+1)
    Q=zeros(T,n,m+1)
    Qp=zeros(T,n,m+1)
    Qn=zeros(T,n,m+1)
    Z=zeros(T,n,m+1)
    H=zeros(T,m+1,m)
    HH=zeros(T,m+1,m)
    ω=zeros(T,m+1)
    a=Vector{T}(γ.^(0:2m+2)); a[1]=zero(T);
    local M0inv::LinSolver = linsolvercreator(nep,σ);
    err=ones(m,m);
    λ=zeros(T,m+1);
    W=zeros(T,n,m+1);

    # temp var for plot
    conv_eig_hist=zeros(Int,m+1)

    # precomputation for exploiting the structure DEP, PEP, GENERAL
    precomp=precompute_data(T,nep,Compute_Bmul_method,n,m,σ,γ)


    # setting initial step
    Q[:,1]=v/norm(v)
    ω[1]=Q[:,1]⋅compute_Mlincomb(nep,0,hcat(Q[:,1],Q[:,1]),[0,1]);
    V[:,1]=Q[:,1];

    k=1; conv_eig=0;
    Av=get_Av(nep)

    # Init errmeasure
    ermdata=init_errmeasure(errmeasure,nep);

    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && ((rem(k,check_error_every)==0) || (k==m))
            println("Iteration:",k, " conveig:",conv_eig)
        end

        broadcast!(/,view(Qn,:,2:k+1),view(Q,:,1:k),(1:k)')
        Qn[:,1] = compute_Mlincomb!(nep,σ,view(Qn,:,1:k+1),a[1:k+1]);
        Qn[:,1] .= -lin_solve(M0inv,Qn[:,1]);

        # call Bmult
        Bmult!(Compute_Bmul_method,k,view(Z,:,1:k+1),Qn,Av,precomp,γ)

        # orthogonalization (three terms recurrence)
        if k>1 β=mat_sum(view(Z,:,1:k),view(Qp,:,1:k)) end
        α=mat_sum(view(Z,:,1:k),view(Q,:,1:k))
        η=mat_sum(view(Z,:,1:k+1),view(Qn,:,1:k+1))

        H[k,k]=α/ω[k]
        if k>1 H[k-1,k]=β/ω[k-1] end
        #Qn[:,1:k] .-= H[k,k]*view(Q,:,1:k);
        mul_and_sub!(view(Qn,:,1:k),view(Q,:,1:k),H[k,k])
        #if k>1 Qn[:,1:k] .-= H[k-1,k]*view(Qp,:,1:k) end
        if k>1 mul_and_sub!(view(Qn,:,1:k),view(Qp,:,1:k),H[k-1,k]) end

        H[k+1,k]=norm(Qn);
#        Qn[:,1:k+1] ./= H[k+1,k]
        scal_mul!(view(Qn,:,1:k+1), 1/H[k+1,k])

        ω[k+1]=η-2*α*H[k,k]+ω[k]*H[k,k]^2;
        if k>1
            ω[k+1]=ω[k+1]-2*β*H[k-1,k]+ω[k-1]*H[k-1,k]^2;
        end
        ω[k+1]=ω[k+1]/H[k+1,k]^2;
        V[:,k+1]=Qn[:,1];

        orthogonalize_and_normalize!(view(V,:,1:k),view(V,:,k+1), view(HH,1:k,k), orthmethod)

        # extract eigenpair approximation
        if ((rem(k,check_error_every)==0)||(k==m))&&(k>2)
            VV=view(V,:,1:k+1)

            # create the projected NEP
            mm=size(VV,2)
            pnep=create_proj_NEP(nep,mm); # maxsize=mm
            set_projectmatrices!(pnep,VV,VV);
            err_lifted=(λ,z)->estimate_error(ermdata,λ,VV*z);

            # solve the projected NEP
            if displaylevel>0
                println("Solving the projected problem")
            end
            λ,ZZ=iar(pnep;Neig=Inf,displaylevel=0,maxit=150,tol=tol,check_error_every=Inf,errmeasure=err_lifted)
            W=VV*ZZ;

            # eigenvectors computation
            conv_eig=length(λ)
            if conv_eig>Neig
                λ=λ[1:Neig]; W=W[:,1:Neig]
            end

        end

        k=k+1;
        # shift the vectors
        Qp=Q;   Q=Qn;
        Qn=zero(Qn);

    end
    k=k-1
    return λ,W,V[:,1:k+1], H[1:k,1:k-1], ω[1:k], HH[1:k,1:k]
end

# this function computes V *= h avoiding allocations (overwrites V)
function scal_mul!(V,h)
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds V[i,j] = h*V[i,j]
        end
    end
end

# this function computes V.-= h*W avoiding allocations (overwrites V)
function mul_and_sub!(V,W,h)
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds V[i,j] -= h*W[i,j]
        end
    end
end

# this function set to zero the variable V avoiding allocations
function mat_zero!(V)
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds V[i,j] = 0
        end
    end
end

# this function computes sum(V.*W) avoiding allocations
function mat_sum(V,W)
    TT=eltype(V)
    β::TT=0
    n,m=size(V)
    for j=1:m
        for i=1:n
            @inbounds β += V[i,j]*W[i,j]
        end
    end
    return β
end

# Contructors for the precomputed data
function PrecomputeDataInit(::Type{Compute_Bmul_method_DEP})
    return IlanPrecomputeDataDEP(0,0,0,0,0)
end
function PrecomputeDataInit(::Type{Compute_Bmul_method_DerSPMF})
    return IlanPrecomputeDataSPMF_NEP(0,0,0,0)
end
function PrecomputeDataInit(::Type{Compute_Bmul_method_SPMF_NEP})
    return IlanPrecomputeDataSPMF_NEP(0,0,0,0)
end

# functions that precompute datas
function precompute_data(T,nep::NEPTypes.DEP,::Type{Compute_Bmul_method_DEP},n,m,σ,γ)
    τ=nep.tauv
    precomp=PrecomputeDataInit(Compute_Bmul_method_DEP)
    precomp.ZZ=zeros(T,n,m+1)       # aux matrix for pre-allocation
    precomp.QQ=zeros(T,n,m+1)       # aux matrix for pre-allocation
    precomp.QQ2=zeros(T,n,m+1)      # aux matrix for pre-allocation
    precomp.vv=zeros(T,m+1,length(τ))
    for j=1:length(τ)
        precomp.vv[:,j]=sqrt(τ[j]*γ)*exp(-σ)*(-τ[j]*γ).^(0:m)
    end
    precomp.G=symmetrizer_coefficients(T,m)

    return precomp
end

# functions that precompute datas
function precompute_data(T,nep::NEPTypes.SPMF_NEP,::Type{Compute_Bmul_method_SPMF_NEP},n,m,σ,γ)
    precomp=PrecomputeDataInit(Compute_Bmul_method_SPMF_NEP)
    precomp.ZZ=zeros(T,n,m+1)       # aux matrix for pre-allocation
    precomp.QQ=zeros(T,n,m+1)       # aux matrix for pre-allocation

    # getting matrices and functions
    fv=get_fv(nep); p=length(fv); Av=get_Av(nep)

    # precompute derivatives and FDH matrices
    SS=diagm(0 => σ*ones(T,2m+2)) + diagm(-1 => γ*(1:2m+1))
    fD=zeros(T,2*m+2,p)
    for t=1:p fD[:,t]=fv[t](SS)[:,1] end
    precomp.FDH=Vector{Array{T,2}}(undef,p)
    for t=1:p precomp.FDH[t]=zeros(T,m+1,m+1)
        for i=1:m+1 for j=1:m+1
            precomp.FDH[t][i,j]=fD[i+j,t];
        end end
    end
    precomp.G=symmetrizer_coefficients(T,m)
    return precomp
end

# functions that precompute datas
function precompute_data(T,nep::NEPTypes.DerSPMF,::Type{Compute_Bmul_method_DerSPMF},n,m,σ,γ)
    precomp=PrecomputeDataInit(Compute_Bmul_method_DerSPMF)
    precomp.ZZ=zeros(T,n,m+1)       # aux matrix for pre-allocation
    precomp.QQ=zeros(T,n,m+1)       # aux matrix for pre-allocation

    # getting matrices and functions
    fv=get_fv(nep); p=length(fv); Av=get_Av(nep)

    precomp.FDH=Vector{Array{T,2}}(undef,p)
    for t=1:p precomp.FDH[t]=zeros(T,m+1,m+1)
        for i=1:m+1 for j=1:m+1
            precomp.FDH[t][i,j]=nep.fD[i+j,t];
        end end
    end
    precomp.G=symmetrizer_coefficients(T,m)
    return precomp
end

function Bmult!(::Type{Compute_Bmul_method_SPMF_NEP},k,Z,Qn,Av,precomp,γ)
    # B-multiplication
    #Z[:,:]=zero(Z);
    mat_zero!(Z)
    @inbounds for t=1:length(Av)
        mul!(view(precomp.QQ,:,1:k+1),view(Qn,:,1:k+1),view(precomp.G,1:k+1,1:k+1).*view(precomp.FDH[t],1:k+1,1:k+1))
        mul!(view(precomp.ZZ,:,1:k+1),Av[t],view(precomp.QQ,:,1:k+1))
        Z .+= view(precomp.ZZ,:,1:k+1);
    end
end

function Bmult!(::Type{Compute_Bmul_method_DerSPMF},k,Z,Qn,Av,precomp,γ)
    Bmult!(Compute_Bmul_method_SPMF_NEP,k,Z,Qn,Av,precomp,γ)
end

function Bmult!(::Type{Compute_Bmul_method_DEP},k,Z,Qn,Av,precomp,γ)
    # B-multiplication
    T=eltype(Z)
    n=size(Z,1)
#    Z[:,:]=zeros(T,n,k+1)
    mat_zero!(Z)
    # low-rank factorization of G
    tolG=1e-12; U,S,V=svd(precomp.G[1:k+1,1:k+1]);q=sum(S.>tolG*ones(length(S)))
    U=view(U,:,1:q).*sqrt.(S[1:q]')
    V=view(V,:,1:q).*sqrt.(S[1:q]');
    Z[:,1]=-γ*Qn[:,1]

    @inbounds for t=2:length(Av)
        mul!(view(precomp.QQ,:,1:q),view(Qn,:,1:k+1),U.*(view(precomp.vv,1:k+1,t-1)))
        # decide multiplication order based on matrix size
        if q>k+1
            mul!(view(precomp.QQ2,:,1:k+1),view(precomp.QQ,:,1:q),(V.*view(precomp.vv,1:k+1,t-1))')
            mul!(view(precomp.ZZ,:,1:k+1),Av[t],view(precomp.QQ2,:,1:k+1));
        else
            mul!(view(precomp.QQ2,:,1:q),Av[t],view(precomp.QQ,:,1:q));
            mul!(view(precomp.ZZ,:,1:k+1),view(precomp.QQ2,:,1:q),(V.*view(precomp.vv,1:k+1,t-1))')
        end
        Z .-= view(precomp.ZZ,:,1:k+1)
    end
end

function symmetrizer_coefficients(T,m)
    G=zeros(T,m+1,m+1); # symmetrizer coefficients
    for i=1:m+1 G[i,1]=1/i end
    for j=1:m for i=1:m+1
            G[i,j+1]=(G[i,j]*j)/(i+j);
    end end
    return G
end
