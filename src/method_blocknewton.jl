using LinearAlgebra

export blocknewton

default_block_errmeasure(nep::NEP) = (S,X) -> opnorm(compute_MM(nep,S,X))

"""
    (S,X)=blocknewton(nep [S,] [X,] [errmeasure,] [tol,] [maxit,] [armijo_factor,] [armijo_max,] [logger])

Applies the block Newton method to `nep::AbstractSPMF`. The method computes
an invariant pair `(S,X)` using the block Newton approach of Kressner.
The variables `S`,`X`
correspond to starting approximations.
The function `errmeasure` shoule be defined for errmeasure(S,X)
and meausures the error in the pair `(S,X)`.
See [`newton`](@ref) for other parameters.


# Example
The example shows that `compute_MM()` becomes zero when a solution
has been computed.
```julia-repl
julia> nep=nep_gallery("dep0",3);
julia> (S,X)= blocknewton(nep)
julia> compute_MM(nep,S,X)
3×2 Array{Complex{Float64},2}:
 -2.22045e-16-1.0842e-19im  -2.08167e-17+0.0im
  1.94289e-16-1.0842e-19im  -5.55112e-17-6.77626e-20im
  7.63278e-17-1.0842e-19im   2.77556e-17-2.71051e-20im
```
This example solves the `gun` problem from the Berlin-Manchester collection
```julia-repl
julia> using NonlinearEigenproblems.Gallery
julia> nep=nep_gallery("nlevp_native_gun");
julia> II=[1.0 0; 0 1]; S=150^2*II; V=[II;zeros(size(nep,1)-2,2)];
julia> (Z,X)=blocknewton(nep,S=S,X=V,logger=1,armijo_factor=0.5,maxit=20)
Iteration 1: Error: 6.081316e+03
Iteration 2: Error: 1.701970e-02 Armijo scaling=0.031250
Iteration 3: Error: 1.814887e-02 Armijo scaling=0.250000
...
Iteration 13: Error: 6.257442e-09
Iteration 14: Error: 2.525942e-15
```

# References
* D. Kressner A block Newton method for nonlinear eigenvalue problems, Numer. Math., 114 (2) (2009), pp. 355-372
"""
function blocknewton(nep::AbstractSPMF;
                     S::Matrix=zeros(2,2),
                     X::Matrix=Matrix(1.0I, size(nep,1), 2),
                     errmeasure::Function =  default_block_errmeasure(nep::NEP),
                     tol::Real=eps(real(eltype(S)))*100,
                     maxit::Integer=10,
                     logger=0,
                     armijo_factor::Real=1,
                     armijo_max::Integer=5)

    @parse_logger_param!(logger)

    T=complex(eltype(S))
    # This implementation is for complex arithmetic only
    # since the original paper is based on the Schur form (not real Schur form)
    # This could potentially be relaxed to real problems if the
    # eigenvalue approximations remain real.
    S=complex(S);
    X=complex(X);

    n=size(nep,1);
    p=size(S,1);

    # Construct the orthogonalization matrix
    W=Vl(X,S);
    # reshape W to WW (3D matrix)
    WW=zeros(T,n,p,p)
    for j=1:size(W,2)
        WW[:,:,j] = W[(j-1)*n .+ (1:n), :]
    end

    local err0=0
    # Main loop
    for k=1:maxit
        err0=errmeasure(S,X)
        push_iteration_info!(logger,k,err=err0,continues=true);
        if (err0<tol)
            push_iteration_info!(logger,"");
            return S,X
        end

        Res= compute_MM(nep,S,X)
        # Solve the linear system by first transforming S to
        # upper triangular form and then doing a backward substitution
        # and finally reversing the backward substitution
        RR,QQ=schur(complex(S))
        dSt,dXt = newtonstep_linsys(T,nep,
                                    RR, X*QQ,
                                    WW,  # Orthogonalization matrix
                                    Res*QQ, # RT
                                    zeros(T,p,p),  #RV=0
                                    logger);
        dX=dXt*QQ';
        dS=QQ*dSt*QQ';


        local j=1
        # Update the approximations
        if (armijo_factor<1)
            (ΔS,ΔV,j,scaling)=  armijo_rule_block(nep,errmeasure,
                                                  err0,
                                                  S,X,-dS,-dX,
                                                  armijo_factor,armijo_max)


            St=S+ΔS
            Xt=X+ΔV
        else
            j = 0
            St=S-dS
            Xt=X-dX

        end
        if (j>0)
            push_info!(logger," Armijo scaling=$scaling")
        else
            push_info!(logger,"");
        end


        # Carry out the orthogonalization
        W,R = qr(Vl(Xt,St))
        W = Matrix(W)
        # reshape W to WW (3D matrix)
        for j=1:size(W,2)
            WW[:,:,j] = W[(j-1)*n .+ (1:n), :]
        end
        X[:]=Xt/R; S[:]=(R*St)/R;
    end
    msg="Number of iterations exceeded. maxit=$(maxit)."
    throw(NoConvergenceException(S,X,err0,msg))

end


# Solve one linear system for a triangular S
# (S,X) Schur pair approx
# W tensor with orthogonalization coeffs
# RT,RV the right-hand side of linear system
function newtonstep_linsys(::Type{T},nep::AbstractSPMF,S, X, W, RT, RV, logger) where {T}

    n = size(nep,1);
    p = size(X,2); l = size(W,3);
    dX = zeros(T, n,p); dS = zeros(T,p,p);

    # Fetch the SPMF
    fv=get_fv(nep);
    Av=get_Av(nep)

    m = size(fv,1);

    # Initialize a tensor with a f_i(S)
    fS = zeros(T,p,p,m);
    for j = 1:m
        fS[:,:,j] = fv[j](S);
    end


    # Work column by-column
    for i = 1:p,
        s = S[i,i];

        # Construct the matrix in equation (20) in Kressner Num. Math.
        T11=compute_Mder(nep,s)
        S_expanded = [S I; zero(S) s*I]
        # T12 = compute_MM(nep,[0*X X],S_expanded) # Can maybe be computed like this? Would avoid the explicit use of Av and and fv
        T12 = zeros(T,n,p);
        for j = 1:m
            DF = fv[j](S_expanded)
            DF1=DF[1:p,p+1:2*p];
            T12 = T12 + Av[j]*X*DF1;
        end
        T21 = W[:,:,1]';
        for j = 2:l
            T21 = T21 + s^(j-1) * W[:,:,j]';
        end
        DS = Matrix(1.0I, p, p)
        T22 = zeros(T,p,p);
        for j = 2:l
            T22 = T22 + W[:,:,j]'*X*DS;
            DS = s*DS + S^(j-2)
        end
        TT=[T11 T12; T21 T22]; # This can maybe be optimized with Schur complement
        sol =  TT \ [RT[:,i];RV[:,i]];
        # compute the dS and dX
        dX[:,i] = sol[1:n];
        dS[:,i] = sol[n+1:end];

        # Update RHS (equation (21)-(22) in Kressner Nummath)
        if (i<p) # Updated RHS needed in the next for loop
            Z = zeros(T,p,p);
            Z[:,i] = dS[:,i]; DS = Z;
            S2_expanded = [S Z; zero(S) S]
            for j = 1:m  # Update\tilde{R}_{T2} according to equation (21)
                Za=dX[:,i]*transpose(fS[i,i+1:p,j]) # First term
                DF = fv[j](S2_expanded)
                Zb=X*DF[1:p,p+i+1:2*p]; # Second term
                RT[:,i+1:p] +=  - Av[j] * (Za+Zb);
            end
            for j = 2:l # Update \tilde{R}_{V2} according to equation (22)
                Za=dX[:,i] * transpose((S^(j-2))[i,i+1:p])
                Zb=X*DS[:,i+1:p] ;
                RV[:,i+1:p] +=  -W[:,:,j]' * ( Za +Zb );
                DS = DS*S + S^(j-2)*DS;
            end
        end
    end
    return dS,dX
end


# Construct the orthogonalization matrix
# V=[X;X*S;XS^2;...]
function Vl(X,S)
    p=size(S,1); n=size(X,1);
    V=zeros(eltype(S),n*p,p);
    for j=1:p
        V[(j-1)*n .+ (1:n), :] = X * (S^(j-1))
    end
    return V
end




function armijo_rule_block(nep,errmeasure,err0,S,V,ΔS,ΔV,armijo_factor,armijo_max)
    j=0
    if (armijo_factor<1)
        # take smaller and smaller steps until errmeasure is decreasing
        while (errmeasure(S+ΔS,V+ΔV)>err0 && j<armijo_max)
            j=j+1;
            ΔS=ΔS*armijo_factor;
            ΔV=ΔV*armijo_factor;
        end
    end
    return  (ΔS,ΔV,j,armijo_factor^j)
end
