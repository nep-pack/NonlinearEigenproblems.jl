using SpecialFunctions
using Random

export infbilanczos
"""
    λv,V,U=infbilanczos([eltype],nep, nept,[linsolvecreator,][linsolvertcreator,][v,][u,][σ,][γ,][tol,][neigs,][errmeasure,][logger,][maxit,][check_error_every])

Executes the Infinite Bi-Lanczos method on the problem defined by `nep::NEP`
and `nept::NEP`. `nep:NEP` is the original nonlinear eigenvalue problem and
`nept::NEP` is its (hermitian) transpose: ``M(λ^*)^H``.
 `v` and `u` are starting vectors,
`σ` is the shift and `γ` the scaling.
The iteration is continued until `neigs` Ritz pairs have converged.
This function throws a `NoConvergenceException` if the wanted eigenpairs are not computed after `maxit` iterations.
However, if `neigs` is set to `Inf` the iteration is continued until `maxit` iterations without an error being thrown.
See [`newton`](@ref) for other parameters.

# Example:
```julia-repl
julia> nep=nep_gallery("dep0");
julia> A=get_Av(nep); fv=get_fv(nep);
julia> At=[copy(A[1]'),copy(A[2]'),copy(A[3]')]
julia> nept=SPMF_NEP(At,fv); # Create the transposed NEP
julia> λv,V=infbilanczos(nep,nept,neigs=3)
julia> norm(compute_Mlincomb(nep,λv[1],V[:,1]))
```

# References:
* The infinite bi-Lanczos method for nonlinear eigenvalue problems, S. W. Gaaf and E. Jarlebring, SIAM J. Sci. Comput. 39:S898-S919, 2017, [preprint](https://arxiv.org/abs/1607.03454)
"""
    infbilanczos(nep::NEP,nept::NEP;params...)=infbilanczos(ComplexF64,nep,nept;params...)
    function infbilanczos(::Type{T},
                          nep::NEP,
                          nept::NEP;  # Transposed NEP
                          maxit::Integer=30,
                          linsolvercreator::Function=default_linsolvercreator,
                          linsolvertcreator::Function=linsolvercreator,
                          v::Vector=randn(real(T),size(nep,1)),
                          u::Vector=randn(real(T),size(nep,1)),
                          tol::Real=1e-12,
                          neigs=5,
                          errmeasure::ErrmeasureType = DefaultErrmeasure,
                          σ::Number=0.0,
                          γ::Number=1,
                          logger=0,
                          check_error_every::Integer=1
                          ) where {T<:Number}

        @parse_logger_param!(logger)

        n=size(nep,1);
        σ=T(σ);
        v=Vector{T}(v);
        u=Vector{T}(v);

        # Linear systems solver for both M(σ) and M(σ)^H

        # Shift σ \neq 0 not implemented

        local M0inv::LinSolver = linsolvercreator(nep,σ);
        local M0Tinv::LinSolver = linsolvertcreator(nept,conj(σ));

        #
        m=maxit;
        qt=lin_solve(M0Tinv,u);
        q=v;
        q=q/dot(qt,compute_Mlincomb(nep,σ,q,ones(1),1));

        Q0=zeros(T,n,m);                  # represents Q_{k-1}
        Qt0=zeros(T,n,m);                 # represents \til{Q}_{k-1}
        R1=zeros(T,n,m+1); R1[:,1]=q;       # represents R_{k}
        Rt1=zeros(T,n,m+1); Rt1[:,1]=qt;    # represents \tild{R}_{k}
        Z2=zeros(T,n,m);
        Zt2=zeros(T,n,m);
        Q_basis=zeros(T,n,m+1);
        Qt_basis=zeros(T,n,m);


        R2=zeros(T,n,m+1); # Needed?
        Rt2=zeros(T,n,m+1); # Needed?

        Q1=zeros(T,n,m); # Needed?
        Qt1=zeros(T,n,m); # Needed?


        # Vectors storing the diagonals
        alpha=zeros(T,m+1);
        beta=zeros(T,m+1);
        gamma=zeros(T,m+1);

        err = zero(T);
        λ=zeros(T,m+1); Q=zeros(T,n,m+1); TT=zeros(T,m+1,m+1);


        normq=norm(q); normqt=norm(qt);
        push_info!(logger,2,"norm(q)=$normq  norm(qt)=$normqt");

        # Init errmeasure
        ermdata=init_errmeasure(errmeasure,nep);


        k=1;

        for k=1:m

            # Note: conjugate required since we compute s'*r not r'*s
            omega = conj(left_right_scalar_prod(T,nep,nept,Rt1,R1,k,k,σ));

            beta[k] = sqrt(abs(omega));
            gamma[k] = conj(omega) / beta[k];

            # Step 11-12

            Q1[:,1:k]=R1[:,1:k]/beta[k];
            Qt1[:,1:k]=Rt1[:,1:k]/conj(gamma[k]);

            #@printf("\n%e %e\n", Q1[1,1],Qt1[1,1]);


            # Extra step, to compute Ritz vectors eventually
            Q_basis[:,k] = Q1[:,1];
            Qt_basis[:,k] = Qt1[:,1];

             # Step 1: Compute Z_{k+1}
            Dk=diagm(0 => 1 ./ (exp.(lfactorial.(1:k))));
            b1_tmp=compute_Mlincomb(nep,σ,Q1[:,1:k]*Dk,ones(k),1);
            b1=-lin_solve(M0inv,b1_tmp);
            Z2[:,k] = b1;



            # Step 2: Compute \til{Z}_{k+1}
            bt1_tmp=compute_Mlincomb(nept,conj(σ),Qt1[:,1:k]*Dk,ones(k),1);
            bt1=-lin_solve(M0Tinv,bt1_tmp);
            Zt2[:,k] = bt1

            # Step 3: Compute R_{k+1}
            R2[:,1] = Z2[:,k];
            R2[:,2:(k+1)]=Q1[:,1:k];
            if k > 1
                R2[:,1:(k-1)]=R2[:,1:(k-1)]-gamma[k]*Q0[:,1:(k-1)];
            end

            # Step 4: Compute \til{R}_{k+1}
            Rt2[:,1] = Zt2[:,k];
            Rt2[:,2:(k+1)]=Qt1[:,1:k];
            if k > 1
                Rt2[:,1:(k-1)]=Rt2[:,1:(k-1)]-conj(beta[k])*Qt0[:,1:(k-1)];
            end

            # Step 5: Compute \alpha_k
            alpha[k+1]=left_right_scalar_prod(T,nep,nept,Qt1,R2,k,k+1,σ);


            #Step 6: Compute R_{k+1}
            R2[:,1:k]=R2[:,1:k]-alpha[k+1]*Q1[:,1:k];



            #Step 7: Compute \til{R}_{k+1}
            Rt2[:,1:k]=Rt2[:,1:k]-conj(alpha[k+1])*Qt1[:,1:k];


            # shift  the matrices:
            (R1,R2)=(R2,R1);  R2[:,:] .= 0  # swap and reset forgotten variable
            (Rt1,Rt2)=(Rt2,Rt1);  Rt2[:,:] .= 0
            (Q0,Q1)=(Q1,Q0);  Q1[:,:] .= 0
            (Qt0,Qt1)=(Qt1,Qt0);  Qt1[:,:] .= 0

            if (rem(k,check_error_every)==0)||(k==m)
                # Check if we should terminate
                omega = left_right_scalar_prod(T,nep,nept,Rt1,R1,k+1,k+1,σ);
                beta[k+1] = sqrt(abs(omega));
                gamma[k+1] = conj(omega) / beta[k+1];

                alpha0=alpha[2:(k+1)];  # \alpha_1 stored in alpha(2)
                beta0=beta[2:(k+1)];    # we do not need \beta_1
                gamma0=gamma[2:(k+1)];  # we do not need \gamma_1

                TT = Matrix(spdiagm(-1 => beta0[1:k], 0 => alpha0[1:k], 1 => gamma0[1:k]))

                E = eigen(TT)
                λ = σ .+ 1 ./ E.values
                Z = E.vectors
                #@ifd(println("size(Z)=",size(Z)))
                #@ifd(println("size(TT)=",size(TT)))
                Q=Q_basis[:,1:(k+1)]*Z
                conv_eig=0;
                # compute the errors

                err=
                  map(s-> estimate_error(ermdata,λ[s],Q[:,s]), 1:size(λ,1))
                # Log them and compute the converged
                push_iteration_info!(logger,2, k,err=err,λ=λ,v=
                                     continues=true);
                for s=1:size(λ,1)
                    if err[s]<tol;
                        conv_eig=conv_eig+1;
                        push_info!(logger,"+", continues=true);
                    elseif err[s]<tol*10
                        push_info!(logger,"=", continues=true);
                    else
                        push_info!(logger,"-", continues=true);
                    end
                end
                push_info!(logger,"");

                idx=sortperm(err[1:k]); # sort the error
                err=err[idx];

                if (conv_eig>=neigs) || (k==m)
                    nrof_eigs = Int(min(length(λ),neigs,conv_eig))
                    λ=λ[idx[1:nrof_eigs]]
                    Q=Q[:,idx[1:nrof_eigs]]
                    for i=1:nrof_eigs
                        normalize!(view(Q,1:n,i))
                    end
                    if (conv_eig>=neigs) || (neigs==Inf)
                        return λ,Q,TT
                    end
                end
            end
        end


        msg="Number of iterations exceeded. maxit=$(maxit)."
        throw(NoConvergenceException(λ,Q,err,msg))
    end

    function left_right_scalar_prod(::Type{T}, nep,nept,At,B,ma,mb,σ) where {T}
        # Compute the scalar product based on the function nep.M_lin_comb
        c=0;
        # This is the nasty double loop, which brings
        # complexity O(m^3n). Will be limiting if we do many iterations
        XX=zeros(T,size(B,1),mb); # pre-allocate
        for j=1:ma
            #dd=1 ./ factorial(j:(j+mb-1));
            dd=1 ./ exp.(lfactorial.(j:(j+mb-1)));
            XX=broadcast(*,B[:,1:mb],dd'); # diag scaling
            #XX=bsxfun(@times,B(:,1:mb),dd);  # Column scaling: Faster than constructing

            # compute Mlincomb starting from derivative j
            z=-compute_Mlincomb(nep,σ,XX,ones(size(XX,2)),j);
            c=c+dot(At[:,j],z);
        end
        return c
    end
