export ilan_benchmark

using IterativeSolvers
using LinearAlgebra
using Random
using Statistics

"""
    ilan_benchmark()
"""
ilan_benchmark(nep::NEP;params...)=ilan_benchmark(ComplexF64,nep;params...)
function ilan_benchmark(
    ::Type{T},
    nep::NEP;
    orthmethod=DGKS(),
    maxit=30,
    linsolvercreator=DefaultLinSolverCreator(),
    tol=eps(real(T))*10000,
    neigs=6,
    errmeasure::ErrmeasureType = DefaultErrmeasure,
    σ=zero(T),
    γ=one(T),
    v=randn(real(T),size(nep,1)),
    displaylevel=0,
    check_error_every=1,
    proj_solve=false,
    inner_solver_method=DefaultInnerSolver())where{T<:Number,T_orth<:IterativeSolvers.OrthogonalizationMethod}

    # Ensure types σ and v are of type T
    σ=T(σ)
    v=Array{T,1}(v)

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
    a=Vector{T}(γ.^(0:m)); a[1]=zero(T); # TODO
    local M0inv::LinSolver = create_linsolver(linsolvercreator,nep,σ)
    err=ones(m,m);
    λ=zeros(T,m+1);

    # precompute the symmetrizer coefficients
    G=zeros(m+1,m+1);
    for i=1:m+1 G[i,1]=1/i end
    for j=1:m
        for i=1:m+1
            G[i,j+1]=(G[i,j]*j)/(i+j);
        end
    end

    # getting matrices and functions
    fv=get_fv(nep); p=length(fv);    Av=get_Av(nep)

    # precompute derivatives and FDH matrices
    SS=diagm(0 => σ*ones(T,2m+2)) + diagm(-1 => (1:2m+1))
    fD=zeros(T,2*m+2,p)
    for t=1:p fD[:,t]=fv[t](SS)[:,1] end
    FDH=Vector{Array{T,2}}(undef,p)
    for t=1:p FDH[t]=zeros(T,m+1,m+1)
        for i=1:m+1 for j=1:m+1
            FDH[t][i,j]=fD[i+j,t];
        end end
    end

    # setting initial step
    Q[:,1]=v/norm(v)
    ω[1]=Q[:,1]⋅compute_Mlincomb(nep,0,hcat(Q[:,1],Q[:,1]),[0,1]);
    V[:,1]=Q[:,1];

    k=1; conv_eig=0;
    local pnep::NEP;
    if (proj_solve)
        pnep=create_proj_NEP(nep);
    end

    while (k <= m) && (conv_eig<neigs)
        if (displaylevel>0) && ((rem(k,check_error_every)==0) || (k==m))
            println("Iteration:",k, " conveig:",conv_eig)
        end

        broadcast!(/,view(Qn,:,2:k+1),view(Q,:,1:k),(1:k)')
        Qn[:,1] = compute_Mlincomb!(nep,σ,Qn[:,1:k+1],a[1:k+1]);
        Qn[:,1] = -lin_solve(M0inv,Qn[:,1]);

        # B-multiplication
        Z=0*Z;
        for t=1:p
            Z[:,1:k+1]=Z[:,1:k+1]+Av[t]*Qn[:,1:k+1]*(G[1:k+1,1:k+1].*FDH[t][1:k+1,1:k+1]);
        end

        # orthogonalization (three terms recurrence)
        α=sum(sum(conj(Z).*Q,dims=1))
        if k>1 β=sum(sum(conj(Z).*Qp,dims=1)) end
        η=sum(sum(conj(Z).*Qn,dims=1))

        H[k,k]=α/ω[k]
        if k>1 H[k-1,k]=β/ω[k-1] end
        Qn=Qn-H[k,k]*Q;
        if k>1 Qn=Qn-H[k-1,k]*Qp end

        H[k+1,k]=norm(Qn);
        Qn=Qn/H[k+1,k];

        ω[k+1]=η-2*α*H[k,k]+ω[k]*H[k,k]^2;
        if k>1
            ω[k+1]=ω[k+1]-2*β*H[k-1,k]+ω[k-1]*H[k-1,k]^2;
        end
        ω[k+1]=ω[k+1]/H[k+1,k]^2;
        V[:,k+1]=Qn[:,1];

        orthogonalize_and_normalize!(view(V,:,1:k),view(V,:,k+1), view(HH,1:k,k), orthmethod)

        k=k+1;
        # shift the vectors
        Qp=Q;   Q=Qn; Qn=0*Qn;

    end
    k=k-1

    return V[:,1:k+1], H[1:k,1:k-1], ω[1:k], HH[1:k,1:k]
end
