export iar

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
iar(nep::NEP;params...)=iar(ComplexF64,nep;params...)
function iar(
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
    V = zeros(T,n*(m+1),m+1);
    H = zeros(T,m+1,m);
    y = zeros(T,n,m+1);
    α=Vector{T}(γ.^(0:m)); α[1]=zero(T);
    local M0inv::LinSolver = linsolvercreator(nep,σ);
    err = NaN*ones(m,m);
    λ=zeros(T,m+1); Q=zeros(T,n,m+1);

    vv=view(V,1:1:n,1); # next vector V[:,k+1]
    vv[:]=v; vv[:]=vv[:]/norm(vv);
    k=1; conv_eig=0;
    local pnep::NEP;
    if (proj_solve)
        pnep=create_proj_NEP(nep);
    end

    # Init errmeasure
    ermdata=init_errmeasure(errmeasure,nep);

    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && ((rem(k,check_error_every)==0) || (k==m))
            println("Iteration:",k, " conveig:",conv_eig)
        end
        VV=view(V,1:1:n*(k+1),1:k); # extact subarrays, memory-CPU efficient
        vv=view(V,1:1:n*(k+1),k+1); # next vector V[:,k+1]

        y[:,2:k+1] = reshape(VV[1:1:n*k,k],n,k);
        broadcast!(/,view(y,:,2:k+1),view(y,:,2:k+1),(1:k)')
        y[:,1] = compute_Mlincomb!(nep,σ,y[:,1:k+1],α[1:k+1]);
        y[:,1] = -lin_solve(M0inv,y[:,1]);

        vv[:]=reshape(y[:,1:k+1],(k+1)*n,1);
        # orthogonalization
        H[k+1,k] = orthogonalize_and_normalize!(VV, vv, view(H,1:k,k), orthmethod)

        # compute Ritz pairs (every check_error_every iterations)
        if ((rem(k,check_error_every)==0)||(k==m))&&(k>2)
            # Extract eigenvalues from Hessenberg matrix
            D,Z = eigen(H[1:k,1:k])

            VV = view(V,1:1:n,1:k)
            Q = VV*Z
            λ = σ .+ γ ./ D

            if (proj_solve)  # Projected solve to extract eigenvalues (otw hessenberg matrix)
                QQ,RR=qr(VV); # Project on this space
                QQ = Matrix(QQ)
                set_projectmatrices!(pnep,QQ,QQ);
                # Make a call to the inner solve method
                λproj,Qproj=inner_solve(inner_solver_method,T,pnep,
                                        V=RR*Z,λv=copy(λ),
                                        Neig=size(λ,1)+3,
                                        σ=mean(λ),
                                        tol=tol,displaylevel=displaylevel);

                Q=QQ*Qproj;
                λ=λproj;
             end

            conv_eig=0;
            for s=1:size(λ,1)
                err[k,s]=estimate_error(ermdata,λ[s],Q[:,s]);
                if err[k,s]<tol; conv_eig=conv_eig+1; end
            end
            idx=sortperm(err[k,1:k]); # sort the error
            err[1:k,k]=err[idx,k];
            # extract the converged Ritzpairs
            if (k==m)||(conv_eig>=Neig)
                if Neig != Inf
                    λ=λ[idx[1:min(length(λ),Neig)]]
                    Q=Q[:,idx[1:length(λ)]]
                else
                    λ=λ[idx[1:min(length(λ))]]
                    Q=Q[:,idx[1:length(λ)]]
                end
            end
        end

        k=k+1;
    end
    k=k-1
    # NoConvergenceException
    if conv_eig<Neig && Neig != Inf
        err=err[end,1:Neig];
        idx=sortperm(err); # sort the error
        λ=λ[idx];  Q=Q[:,idx]; err=err[idx];
        msg="Number of iterations exceeded. maxit=$(maxit)."
        if conv_eig<3
            msg=string(msg, " Check that σ is not an eigenvalue.")
        end
        throw(NoConvergenceException(λ,Q,err,msg))
    end

    # extract the converged Ritzpairs
    λ=λ[1:min(length(λ),conv_eig)];
    Q=Q[:,1:min(size(Q,2),conv_eig)];
    return λ,Q,err[1:k,:],V[:,1:k]
end
