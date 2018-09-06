export ilan

using IterativeSolvers
using Random

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
ilan(nep::NEP;params...)=ilan(ComplexF64,nep;params...)
function ilan(
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
    inner_solver_method=DefaultInnerSolver) where {T,T_orth<:IterativeSolvers.OrthogonalizationMethod}


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

    k=1; conv_eig=0;

    V[1:n,2]=v;
    ω=zeros(T,m+1); # normalization coefficients
    ω[1]=V[:,2]⋅B_action(T,nep,V[:,2]);

    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && ((rem(k,check_error_every)==0) || (k==m))
            println("Iteration:",k, " conveig:",conv_eig)
        end

        # action of C in the last vector (identical to iar)
        y[:,2:k+1] = reshape(V[1:1:n*k,2],n,k);
        broadcast!(/,view(y,:,2:k+1),view(y,:,2:k+1),(1:k)')
        y[:,1] = compute_Mlincomb(nep,σ,y[:,1:k+1],α[1:k+1]);
        y[:,1] = -lin_solve(M0inv,y[:,1]);
        V[1:(k+1)*n,3]=reshape(y[:,1:k+1],(k+1)*n,1);

        # orthogonalization
        z=B_action(T,nep,V[:,3]);
        alpha=z⋅V[:,2]; gamma=z⋅V[:,3];
        H[k,k]=alpha/ω[k];
        V[:,3]=V[:,3]-H[k,k]*V[:,2];
        ω[k+1]=gamma-2*H[k,k]*alpha+H[k,k]^2*ω[k];

        if k>1
            beta=z⋅V[:,1];
            H[k-1,k]=beta/ω[k-1];
            V[:,3]=V[:,3]-H[k-1,k]*V[:,1];
            ω[k+1]=ω[k+1]-2*H[k-1,k]*beta+H[k-1,k]^2*ω[k-1];
        end

        # shift the three term recurrence
        V[:,1:2]=V[:,2:3];

        k=k+1;
    end
    k=k-1


    return V
end

function B_action(T,nep::NEPTypes.NEP,z)
    return z
end
