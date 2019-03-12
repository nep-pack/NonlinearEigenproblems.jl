export ilan_benchmark

using IterativeSolvers
using LinearAlgebra
using Random
using Statistics

"""
    iar(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=1,][orthmethod=DGKS])

Run the infinite Arnoldi method on the nonlinear eigenvalue problem stored in `nep`.

The target `σ` is the center around which eiganvalues are computed. The kwarg `errmeasure` is a function handle which can be used to specify how the error is measured to be used in termination (default is absolute residual norm). A Ritz pair `λ` and `v` is flagged a as converged (to an eigenpair) if `errmeasure` is less than `tol`. The vector
`v` is the starting vector for constructing the Krylov space. The orthogonalization method, used in contructing the orthogonal basis of the Krylov space, is specified by `orthmethod`, see the package `IterativeSolvers.jl`. The iteration
is continued until `Neig` Ritz pairs converge. This function throws a `NoConvergenceException` if the wanted eigenpairs are not computed after `maxit` iterations. The `linsolvercreator` is a function which specifies how the linear system is created and solved.


# Example
```julia-repl
julia> using NonlinearEigenproblems, LinearAlgebra
julia> nep=nep_gallery("dep0",100);
julia> v0=ones(size(nep,1));
julia> λ,v=iar(nep;v=v0,tol=1e-5,Neig=3);
julia> norm(compute_Mlincomb!(nep,λ[1],v[:,1])) # Is it an eigenvalue?
julia> λ    # print the computed eigenvalues
3-element Array{Complex{Float64},1}:
 -0.15606211475666945 - 0.12273439802763578im
 -0.15606211475666862 + 0.12273439802763489im
  0.23169243065648365 - 9.464790582509696e-17im
```

# References
* Algorithm 2 in Jarlebring, Michiels Meerbergen, A linear eigenvalue algorithm for the nonlinear eigenvalue problem, Numer. Math, 2012
"""
ilan_benchmark(nep::NEP;params...)=ilan_benchmark(ComplexF64,nep;params...)
function ilan_benchmark(
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
    check_error_every=1,
    proj_solve=false,
    inner_solver_method=DefaultInnerSolver)where{T<:Number,T_orth<:IterativeSolvers.OrthogonalizationMethod}

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
    local M0inv::LinSolver = linsolvercreator(nep,σ);
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

    while (k <= m) && (conv_eig<Neig)
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
