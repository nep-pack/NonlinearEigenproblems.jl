export ilan
using IterativeSolvers

"""
    ilan(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure=default_errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=1,][orthmethod=DGKS])

### Infinite Arnoldi method

Runs the infinite Arnoldi method which tries to find eigenvalues close to the shift σ.


# Example
```julia-repl
julia> using NonlinearEigenproblems: NEPSolver, NEPCore, Gallery
julia> nep=nep_gallery("dep0");
julia> λ,v=iar(nep);
julia> minimum(svdvals(compute_Mder(nep,λ[1]))) % Is it an eigenvalue?

```

# References
* Algorithm 2 in Jarlebring, Michiels Meerbergen, A linear eigenvalue algorithm for the nonlinear eigenvalue problem, Numer. Math, 2012
"""
ilan(nep::NEP;params...)=ilan(Complex128,nep;params...)
function ilan{T,T_orth<:IterativeSolvers.OrthogonalizationMethod}(
    ::Type{T},
    nep::NEP;
    orthmethod::Type{T_orth}=DGKS,
    maxit=30,
    linsolvercreator::Function=default_linsolvercreator,
    tol=eps(real(T))*10000,
    Neig=6,
    errmeasure::Function = default_errmeasure(nep::NEP),
    σ=zero(T),
    γ=one(T),
    v=randn(real(T),size(nep,1)),
    displaylevel=0,
    check_error_every=1,
    proj_solve=false,
    inner_solver_method=DefaultInnerSolver)


    n = size(nep,1);
    m = maxit;

    # initialization

    V = zeros(T,n*(m+1),3); # three term recurrence
    H = zeros(T,m+1,m);     # projection matrix
    y = zeros(T,n,m+1);     # auxiliary vector

    α=γ.^(0:m); α[1]=zero(T);
    local M0inv::LinSolver = linsolvercreator(nep,σ);
    err = ones(m,m);
    λ=zeros(T,m+1); Q=zeros(T,n,m+1);

    vv=view(V,1:1:n,1); # next vector V[:,k+1]
    vv[:]=v; vv[:]=vv[:]/norm(vv);
    k=1; conv_eig=0;
    local pnep::NEP;
    if (proj_solve)
        pnep=create_proj_NEP(nep);
    end

    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && ((rem(k,check_error_every)==0) || (k==m))
            println("Iteration:",k, " conveig:",conv_eig)
        end

        # action of C in the last vector (identical to iar)
        y[:,2:k+1] = reshape(V[1:1:n*k,2],n,k);
        broadcast!(/,view(y,:,2:k+1),view(y,:,2:k+1),(1:k)')
        y[:,1] = compute_Mlincomb(nep,σ,y[:,1:k+1],a=α[1:k+1]);
        y[:,1] = -lin_solve(M0inv,y[:,1]);

        V[1:(k+1)*n,3]=reshape(y[:,1:k+1],(k+1)*n,1);

        # orthogonalization
        # missing



        k=k+1;
    end
    k=k-1


    return V
end
