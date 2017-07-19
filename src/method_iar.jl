export iar
# The Infinite Arnoldi method

"""
### Infinite Arnoldi method
Compute `Neig` eigenvalues closer to `σ` and the corresponing eigenvectors. 
"""

function iar(
    nep::NEP;
    maxit=30,
    linsolvercreator::Function=default_linsolvercreator,
    tol=1e-12,
    Neig=maxit,
    errmeasure::Function = default_errmeasure(nep::NEP),
    σ=0.0,
    γ=1,
    v0=rand(size(nep,1),1),
    displaylevel=0,
    p=1
    )

    n = size(nep,1);
    m = maxit;

    # initialization
    V = zeros(Complex128,n*(m+1),m+1);
    H = zeros(Complex128,m+1,m);
    y = zeros(Complex128,n,m+1);
    α=γ.^(0:m); α[1]=0;
    local M0inv::LinSolver = linsolvercreator(nep,σ);
    err = ones(m,m);
    λ=complex(zeros(m+1)); Q=complex(zeros(n,m+1));

    vv=view(V,1:1:n,1); # next vector V[:,k+1]
    vv[:]=v0; vv[:]=vv[:]/norm(vv);
    k=1; conv_eig=0;

    # check if σ is a good approximation of an eigenvalue
    v=compute_eigvec_from_eigval(nep,λ[1], (nep, σ) -> M0inv);
    if errmeasure(σ,v)<tol
        msg="The shift σ is an eigenvalue, try a different value for σ."
        throw(NoConvergenceException(σ,v,err,msg))
    end

    while (k <= m) && (conv_eig<Neig)
        if (displaylevel>0) && (rem(k,p)==0) || (k==m)
         println("Iteration:",k, " conveig:",conv_eig)
    end
    VV=view(V,1:1:n*(k+1),1:k); # extact subarrays, memory-CPU efficient
    vv=view(V,1:1:n*(k+1),k+1); # next vector V[:,k+1]

    y[:,2:k+1] = reshape(VV[1:1:n*k,k],n,k);
    for j=1:k
        y[:,j+1]=y[:,j+1]/j;
    end
    y[:,1] = compute_Mlincomb(nep,σ,y[:,1:k+1],a=α[1:k+1]);
    y[:,1] = -lin_solve(M0inv,y[:,1]);

    vv[:]=reshape(y[:,1:k+1],(k+1)*n,1);
    # orthogonalization
    h,vv[:] = doubleGS(VV,vv,k,n);
    H[1:k,k]=h;
    β=norm(vv);

    H[k+1,k]=β;
    vv[:]=vv[:]/β;

    # compute Ritz pairs (every p iterations)
    if (rem(k,p)==0)||(k==m)
        D,Z=eig(H[1:k,1:k]);
        VV=view(V,1:1:n,1:k);
        Q=VV*Z; λ=σ+γ./D;
        conv_eig=0;
        for s=1:k
            err[k,s]=errmeasure(λ[s],Q[:,s]);
            if err[k,s]>10; err[k,s]=1; end	# artificial fix
            if err[k,s]<tol; conv_eig=conv_eig+1; end
        end
        idx=sortperm(err[1:k,k]); # sort the error
        err[1:k,k]=err[idx,k];
        # extract the converged Ritzpairs
        if (k==m)||(conv_eig>=Neig)
            λ=λ[idx[1:min(length(λ),Neig)]]
            Q=Q[:,idx[1:length(λ)]]
        end
        end
    k=k+1;
    end

    if conv_eig<Neig
        err=err[end,1:Neig];
        idx=sortperm(err); # sort the error
        λ=λ[idx];  Q=Q[:,idx]; err=err[idx];
        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,Q,err,msg))
    end

    return λ,Q,err
end


function doubleGS(VV,vv,k,n)
    h=VV'*vv;
    vv=vv-VV*h;
    g=VV'*vv;
    vv=vv-VV*g;
    h = h+g;
    return h,vv;
end

function singleGS(V,vv,k,n)
    h=V[1:(k+1)*n,1:k]'*vv;
    vv=vv-V[1:(k+1)*n,1:k]*h;
    return h,vv;
end
