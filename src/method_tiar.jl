export tiar

using IterativeSolvers
using LinearAlgebra
using Random

"""
    tiar(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=1,][orthmethod=DGKS])

Run the tensor infinite Arnoldi method on the nonlinear eigenvalue problem stored in `nep`.

The target `σ` is the center around which eiganvalues are computed. The kwarg `errmeasure` is a function handle which can be used to specify how the error is measured to be used in termination (default is absolute residual norm). A Ritz pair `λ` and `v` is flagged a as converged (to an eigenpair) if `errmeasure` is less than `tol`. The vector
`v` is the starting vector for constructing the Krylov space. The orthogonalization method, used in contructing the orthogonal basis of the Krylov space, is specified by `orthmethod`, see the package `IterativeSolvers.jl`. The iteration
is continued until `Neig` Ritz pairs converge. This function throws a `NoConvergenceException` if the wanted eigenpairs are not computed after `maxit` iterations. The `linsolvercreator` is a function which specifies how the linear system is created and solved.

# Example
```julia-repl
julia> using NonlinearEigenproblems, LinearAlgebra
julia> nep=nep_gallery("dep0",100);
julia> v0=ones(size(nep,1));
julia> λ,v=tiar(nep;v=v0,tol=1e-5,Neig=3);
julia> norm(compute_Mlincomb!(nep,λ[1],v[:,1])) # Is it an eigenvalue?
julia> λ    # print the computed eigenvalues
3-element Array{Complex{Float64},1}:
 -0.1560621147566685 + 0.12273439802763504im
 -0.1560621147566693 - 0.1227343980276357im
 0.23169243065648332 - 4.699260229885766e-17im

```

# References
* Algorithm 2 in Jarlebring, Mele, Runborg, The Waveguide Eigenvalue Problem and the Tensor Infinite Arnoldi Method, SIAM J. Scient. computing, 39 (3), A1062-A1088, 2017

"""
tiar(nep::NEP;params...)=tiar(ComplexF64,nep;params...)
function tiar(
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
    inner_solver_method=DefaultInnerSolver
    )where{T,T_orth<:IterativeSolvers.OrthogonalizationMethod}

    # Ensure types σ and v are of type T
    σ=T(σ)
    v=Array{T,1}(v)

    # initialization
    n = size(nep,1); m = maxit;

    if n<m
        msg="Loss of orthogonality in the matrix Z. The problem size is too small, use iar instead.";
        throw(LostOrthogonalityException(msg))
    end

    # initialize variables
    a  = zeros(T,m+1,m+1,m+1);
    Z  = zeros(T,n,m+1);
    t  = zeros(T,m+1);
    tt = zeros(T,m+1);
    g  = zeros(T,m+1,m+1);
    f  = zeros(T,m+1,m+1);
    ff = zeros(T,m+1,m+1);
    H  = zeros(T,m+1,m);
    h  = zeros(T,m+1);
    hh = zeros(T,m+1);
    y  = zeros(T,n,m+1);
    α=Array{T,1}(γ.^(0:m)); α[1]=zero(T);
    local M0inv::LinSolver = linsolvercreator(nep,σ);
    err = NaN*ones(m+1,m+1);
    λ=zeros(T,m+1); Q=zeros(T,n,m+1);
    Z[:,1]=v; Z[:,1]=Z[:,1]/norm(Z[:,1]);
    a[1,1,1]=one(T);

    # temp var for plot
    conv_eig_hist=zeros(Int,m+1)

    local pnep::NEP;
    if (proj_solve)
        pnep=create_proj_NEP(nep,maxit,T);
    end

    # Init errmeasure
    ermdata=init_errmeasure(errmeasure,nep);

    k=1; conv_eig=0;
    while (k <= m)&(conv_eig<Neig)
        if (displaylevel>0) && ((rem(k,check_error_every)==0) || (k==m))
            println("Iteration:",k, " conveig:",conv_eig)
        end

        # computation of y[:,2], ..., y[:,k+1]
        y[:,2:k+1]=Z[:,1:k]*transpose(a[1:k,k,1:k])
        broadcast!(/,view(y,:,2:k+1),view(y,:,2:k+1),(1:k)')

        # computation of y[:,1]
        y[:,1] = compute_Mlincomb!(nep,σ,y[:,1:k+1],α[1:k+1]);
        y[:,1] = -lin_solve(M0inv,y[:,1]);

        # Gram–Schmidt orthogonalization in Z
        Z[:,k+1]=y[:,1];
        t[k+1] = orthogonalize_and_normalize!(view(Z,:,1:k), view(Z,:,k+1), view(t,1:k), orthmethod)

        # compute the matrix G
        for l=1:k+1
            for i=2:k+1
                g[i,l]=a[i-1,k,l]/(i-1);
            end
         g[1,l]=t[l];
        end

        # compute h (orthogonalization with tensors factorization)
        h = zero(h)
        Ag = zero(h[1:k])
        for l=1:k
            mul!(Ag,a[1:k,1:k,l]',g[1:k,l])
            h[1:k] .+= Ag;
        end

        # compute the matrix F
        f=g;
        Ah = zero(f[1:k+1,1])
        for l=1:k
            mul!(Ah,a[1:k+1,1:k,l],h[1:k])
            f[1:k+1,l] .-= Ah;
        end

        # re-orthogonalization
        # compute hh (re-orthogonalization with tensors factorization)
        hh = zero(hh)
        Af = zero(hh[1:k])
        for l=1:k
            mul!(Af,a[1:k,1:k,l]',f[1:k,l])
            hh[1:k] .+= Af;
        end

        # compute the matrix FF
        ff=f;
        Ah=zero(ff[1:k+1,1])
        for l=1:k
            mul!(Ah,a[1:k+1,1:k,l],hh[1:k])
            ff[1:k+1,l] .-= Ah;
        end

        # update the orthogonalization coefficients
        h=h+hh; f=ff;
        β=norm(view(f,1:k+1,1:k+1)); # equivalent to Frobenius norm

        # extend the matrix H
        H[1:k,k]=h[1:k]; H[k+1,k]=β;

        # extend the tensor
        for i=1:k+1
            for l=1:k+1
                a[i,k+1,l]=f[i,l]/β;
            end
        end

        # compute Ritz pairs (every p iterations)
        if (rem(k,check_error_every)==0)||(k==m)
            D,W = eigen(H[1:k,1:k])
            VV=Z[:,1:k]*transpose(a[1,1:k,1:k])	# extract proper subarray
            Q = VV*W
            λ = σ .+ γ ./ D

            if (proj_solve)  # Projected solve to extract eigenvalues (otw hessenberg matrix)
                set_projectmatrices!(pnep,Z[:,1:k],Z[:,1:k]);
                # Make a call to the inner solve method
                λproj,Qproj=inner_solve(inner_solver_method,T,pnep,
                                        λv=copy(λ),
                                        Neig=size(λ,1)+3,
                                        σ=σ,
                                        tol=tol/10,displaylevel=displaylevel);


                II=sortperm(abs.(λproj .- σ));
                λproj=λproj[II]; Qproj=Qproj[:,II];
                Q=Z[:,1:k]*Qproj;
                λ=λproj;
             end


            conv_eig=0;
            for s=1:min(size(λ,1),size(err,2))
                err[k,s]=estimate_error(ermdata,λ[s],Q[:,s]);
                if err[k,s]<tol; conv_eig=conv_eig+1; end
            end
            idx=sortperm(err[k,1:k]); # sort the error
            err[1:k,k]=err[idx,k];
            # extract the converged Ritzpairs
            if (k==m)||(conv_eig>=Neig)
                λ=λ[idx[1:min(length(λ),Neig)]]
                Q=Q[:,idx[1:length(λ)]]
            end
            conv_eig_hist[k]=conv_eig
        end
        k=k+1;
    end

    # NoConvergenceException
    # if conv_eig<Neig
    #    err=err[end,1:Neig];
    #    idx=sortperm(err); # sort the error
    #    λ=λ[idx];  Q=Q[:,idx]; err=err[idx];
    #     msg="Number of iterations exceeded. maxit=$(maxit)."
    #     if conv_eig<3
    #         msg=string(msg, " Check that σ is not an eigenvalue.")
    #     end
    #     throw(NoConvergenceException(λ,Q,err,msg))
    # end

    # extract the converged Ritzpairs
    λ=λ[1:min(length(λ),conv_eig)];
    Q=Q[:,1:min(size(Q,2),conv_eig)];
    return λ,Q,err[1:k,:],Z[:,1:k],conv_eig_hist
end
