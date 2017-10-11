export iar_chebyshev
using IterativeSolvers


"""
    iar(nep,[maxit=30,][σ=0,][γ=1,][linsolvecreator=default_linsolvecreator,][tolerance=eps()*10000,][Neig=6,][errmeasure=default_errmeasure,][v=rand(size(nep,1),1),][displaylevel=0,][check_error_every=1,][orthmethod=DGKS])
### Infinite Arnoldi method
Infinite Arnoldi method, as described in Algorithm 2 in  "A linear eigenvalue algorithm for the nonlinear eigenvalue problem",
by Jarlebring, Elias and Michiels, Wim and Meerbergen, Karl.
"""
iar_chebyshev(nep::NEP;params...)=iar_chebyshev(Complex128,nep;params...)
function iar_chebyshev{T,T_orth<:IterativeSolvers.OrthogonalizationMethod}(
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
    check_error_every=1
    )

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
    println("Type in iar", typeof(λ))

    vv=view(V,1:1:n,1); # next vector V[:,k+1]
    vv[:]=v; vv[:]=vv[:]/norm(vv);
    k=1; conv_eig=0;

    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && (rem(k,check_error_every)==0) || (k==m)
            println("Iteration:",k, " conveig:",conv_eig)
        end
        VV=view(V,1:1:n*(k+1),1:k); # extact subarrays, memory-CPU efficient
        vv=view(V,1:1:n*(k+1),k+1); # next vector V[:,k+1]

        # just a test
        # hardcoded matrix L
        temp=reshape(vv[1:k*n],n,k);
        println("temp",temp)

        if k==1
            L=2;
        else
            L=diagm(vcat(2, 1./(2:k)),0)+diagm(-vcat(1./(1:(k-2))),-2);
        end
        y[:,2:k+1]=reshape(vv[1:k*n],n,k)*L;
        # end test
        #println("Matrix y",y)

        y[:,1] = compute_y0(reshape(vv[1:k*n],n,k),y[:,2:k+1],nep);
        #println("Matrix y after y0",y)

        vv[:]=reshape(y[:,1:k+1],(k+1)*n,1);
        # orthogonalization
        H[k+1,k] = orthogonalize_and_normalize!(VV, vv, view(H,1:k,k), orthmethod)

        #println("Vector vv",vv)
        #println("Matrix H",H[1:k,1:k])

        # compute Ritz pairs (every check_error_every iterations)
        if (rem(k,check_error_every)==0)||(k==m)
            D,Z=eig(H[1:k,1:k]);
            VV=view(V,1:1:n,1:k);
            Q=VV*Z; λ=D;
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

    return λ,Q,err[1:k,:],V[:,1:k]
end


function compute_y0(x,y,nep)
   n,N=size(y);
   y0=zeros(n,1);
   A0=nep.A[2];
   A1=nep.A[3];


   for i=1:N
      y0=y0+T(i,3)*y[:,i];
   end
   y0=-A1*y0;

   for i=1:N-1
      # be careful with the index of x. It starts from zero.
      y0=y0+2*i*U(i-1,1)*x[:,i+1];
   end
   y0=y0-A0*sum(y,2);
   y0=(A0+A1)\y0;

end


# Chebyshev polynomials of the first kind
function T(i,x)
   if i==0
      1
   elseif i==1
      x
   else
      2*x*T(i-1,x)-T(i-2,x)
   end
end

# Chebyshev polynomials of the second kind
function U(i,x)
   if i==0
      1
   elseif i==1
      2*x
   else
      2*x*U(i-1,x)-U(i-2,x)
   end
end
