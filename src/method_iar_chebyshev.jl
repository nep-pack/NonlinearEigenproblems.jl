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
    errmeasure::Function=default_errmeasure(nep::NEP),
    σ=zero(T),
    γ=one(T),
    v=randn(real(T),size(nep,1)),
    displaylevel=0,
    check_error_every=1,
    compute_y0::Function=function emptyfunc end
    )
    # hardcoded for 2dep
    a=-5; b=3;

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

    vv=view(V,1:1:n,1); # next vector V[:,k+1]
    v=ones(n,1);        # debug
    vv[:]=v; vv[:]=vv[:]/norm(vv);
    k=1; conv_eig=0;

    # hardcoded matrix L
    L=diagm(vcat(2, 1./(2:m)),0)+diagm(-vcat(1./(1:(m-2))),-2);
    L=L*(b-a)/4;

    # Compute the P and P_inv
    P=P_mat(T,m+1,2/(b-a),(a+b)/(a-b));
    P_inv=P_inv_mat(T,m+1,2/(b-a),(a+b)/(a-b));
    #println("\n")
    #Base.showarray(STDOUT,P,false)
    #println("\n")
    #Base.showarray(STDOUT,P_inv,false)
    #println("\n")


    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && (rem(k,check_error_every)==0) || (k==m)
            println("Iteration:",k, " conveig:",conv_eig)
        end
        VV=view(V,1:1:n*(k+1),1:k); # extact subarrays, memory-CPU efficient
        vv=view(V,1:1:n*(k+1),k+1); # next vector V[:,k+1]

        # compute y (only steps different then standard IAR)
        y=zeros(n,k+1);
        if isa(nep,NEPTypes.DEP)
            y[:,2:k+1]  = reshape(VV[1:1:n*k,k],n,k)*L[1:k,1:k];
            y[:,1]      = compute_y0_dep(reshape(VV[1:1:n*k,k],n,k),y[:,1:k+1],nep,a,b);
        elseif isempty(methods(compute_y0))
            y[:,2:k+1] = reshape(VV[1:1:n*k,k],n,k)*P[1:k,1:k]';
            for j=1:k
                y[:,j+1]=y[:,j+1]/j;
            end
            y[:,1] = compute_Mlincomb(nep,σ,y[:,1:k+1],a=α[1:k+1]);
            y[:,1] = -lin_solve(M0inv,y[:,1]);

            y=y*P_inv[1:k+1,1:k+1]';
        else
            y[:,2:k+1]  = reshape(VV[1:1:n*k,k],n,k)*L[1:k,1:k];
            y[:,1]      = compute_y0(reshape(VV[1:1:n*k,k],n,k),y[:,1:k+1],nep,a,b);
        end

        vv[:]=reshape(y[:,1:k+1],(k+1)*n,1);
        # orthogonalization
        H[k+1,k] = orthogonalize_and_normalize!(VV, vv, view(H,1:k,k), orthmethod)

        # compute Ritz pairs (every check_error_every iterations)
        if (rem(k,check_error_every)==0)||(k==m)
            D,Z=eig(H[1:k,1:k]);
            VV=view(V,1:1:n,1:k);
            Q=VV*Z; λ=1./D;
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

    k=k-1
    return λ,Q,err[1:k,:],V,H
end



function mon2cheb(T,ρ,γ,a)
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


function P_mat(T,n, ρ, γ )
    I=eye(T,n,n); P=zeros(T,n,n);
    ρ=T(ρ); γ=T(γ);
    for j=1:n
        P[:,j]=cheb2mon(T,ρ,γ,I[:,j]);
    end
    return P
end


function P_inv_mat(T,n, ρ, γ )
    I=eye(T,n,n); P_inv=zeros(T,n,n);
    ρ=T(ρ); γ=T(γ);
    for j=1:n
        P_inv[:,j]=mon2cheb(T,ρ,γ,I[:,j]);
    end
    return P_inv
end


function compute_y0_dep(x,y,nep,a,b)
   # TODO: DOCUMENT THIS BETTER
   #y_0= \sum_{i=1}^N T_{i-1}(c) x_i - \sum_{j=1}^m A_j \left( \sum_{i=1}^{N+1} T_{i-1}(-k \tau_j+c) y_i\right )
   T=(n,x)->cos(n*acos(x));
   c=(a+b)/(a-b);   k=2/(b-a);

   N=size(x,2);
   n=size(x,1);


   y0=zeros(n,1);
   # TODO: IMPROVE THESE FOR LOOPS
   for i=1:N
       y0=y0+T(i-1,c)*x[:,i];
   end

   # TODO: IMPROVE THESE FOR LOOPS
   for j=1:length(nep.A)
       for i=1:N+1
           y0=y0-T(i-1,-k*nep.tauv[j]+c)*(nep.A[j]*y[:,i]);
       end
   end

   # USE THE LINSOLVER
   y0=sum(nep.A)\y0;
   return y0
end
