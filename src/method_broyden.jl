using LinearAlgebra
#using Arpack # disabled. Use ArnoldMethod instead
using ..NEPCore, ..NEPTypes
import ..NEPCore.compute_Mder;
import ..NEPCore.compute_Mlincomb;
import Base.size;
export broyden;





function broyden_default_errmeasure(λ,v,r)
    return norm(r)/norm(v)
end



# The internal Broyden method (without deflation)
# which corresponds to version T in the paper below
function broyden_T(::Type{TT},nep::NEP;
                   v1=0,u1=[],λ1=TT(0),
                   CH=0,T1=0,W1=0,
                   maxit=100,
                   S=zeros(TT,0,0),X=zeros(TT,size(nep,1),0),
                   check_error_every=10,
                   print_error_every=1,
                   tol=1e-12,
                   threshold=0.4,
                   time0=time_ns(),
                   errmeasure::Function=broyden_default_errmeasure,
                   logger=0) where {TT<:Number}

    @parse_logger_param!(logger)

    # Check types are consistent

    v=Array{TT,1}(v1);
    u=Array{TT,1}(u1);
    X=Array{TT,2}(X);
    S=Array{TT,2}(S);
    CH=Array{TT,2}(CH);
    λ=TT(λ1);
    #x=vcat(v,λ);
    #x=Array{TT,1}(x);
    n=size(nep,1);
    p=size(S,1);

    rk = compute_Mlincomb(nep, λ, v+(X*((λ*Matrix{TT}(I, p, p) - S) \ u)))

    T = Matrix{TT}(T1)
    Wext = Matrix{TT}(undef, n, p+2)
    W = view(Wext, 1:n, 1:p+1)
    W[:,:] = Matrix{TT}(W1)
    errhist = Vector{real(TT)}(NaN * ones(real(TT), maxit))

    timehist = Vector{Float64}(undef, maxit) * NaN
    II = Matrix{TT}(I, p, p)

    Z = Matrix{TT}(undef, n, p+1)
    Tztilde = Vector{TT}(undef, n)
    aH = Vector{TT}(undef, p)'
    Wold = copy(W)
    ztilde = Vector{TT}(undef, n)
    waH = Vector{TT}(undef, p)'
    Z=T*W;
    for j=1:maxit

        Trk=T*rk;

        duλ=-(CH*Z)\(CH*Trk) # contains [Δu;Δλ]


        # Compute u and λ updates
        Δu=duλ[1:p];
        Δλ=duλ[end];
        # Compute v-update
        Δv=-Z*duλ-Trk



        γ=TT(1.0);
        tt=real(TT)(sqrt(abs(Δλ)^2+norm(Δv)^2));
        if (tt>threshold)  # Avoid too big λ-steps
            γ=TT(threshold)/tt;
        end


        v=v+γ*Δv;
        u=u+γ*Δu;
        λ=λ+γ*Δλ;



        # Update Jacobian approximation
        vv=v+(X*((λ*II-S)\u));
        rkp::AbstractVector=compute_Mlincomb(nep,λ,vv);


        ztilde=(rkp-(1-γ)*rk)/γ;
        Tztilde=T*ztilde;

        # For w-update
        bH=[Δu' (Δλ)']/(norm(Δv)^2+norm(Δu)^2+norm(Δλ)^2);
        # For T-update
        β=norm(Δv)^2+norm(Δu)^2+abs(Δλ)^2+Δv'*Tztilde;
        aH=-(Δv'*T)/β

        # Update Z
        Z=Z+Tztilde*(aH*W+(1+aH*ztilde)*bH);

        # Update W
        W=W+ztilde*bH;


        # Update T
        T=T+((Tztilde)*aH);

        if (false) # Extra type check
            for q=[Z,T,W,Trk,duλ,Δu,Δv,λ,Δλ,Tztilde,T,aH,complex(β)]
                if (TT != eltype(q))
                    error("Bad type ",q);
                end

            end
        end





        # step forward
        rk=rkp;

        if (mod(j,check_error_every)==0)
            errhist[j]=errmeasure(λ,v+(X/(λ*II-S))*u,rk);
            timehist[j]=Float64((time_ns()-time0)*1e-9);
            if (mod(j,print_error_every)==0)
                d = opnorm(CH*v - reverse(Matrix{TT}(I, 1+p, 1), dims = 1))
            end



            push_iteration_info!(logger,j,err=errhist[j],λ=λ);
            if (errhist[j]<tol)
                return (λ,v,u,T,W,j,errhist[1:j],timehist[1:j])
            end
        end

    end

    push_info!(logger,"Too many iterations"); #" resnorm=",opnorm(rk));
    return (λ,v,u,T,W,maxit,errhist,timehist)

end

##################
"""
    S,V = broyden([eltype,]nep::NEP[,approxnep];kwargs)

Runs Broydens method (with deflation) for the nonlinear eigenvalue
problem defined by nep.
An approximate nep can be provided which is used as an
initialization of starting
matrix/vectors. The optional argument `approxnep` determines how
 to initiate the
algorithm. It can be an `NEP`, the symbol `:eye`
corresponding to starting
with an identity matrix, and a `Matrix` (of size ``n\times n``).
Beside most of the standard kwargs as described in [`augnewton`](@ref),
it supports `pmax` which is subspace used in deflation, essentially
the number of eigenvalues, `add_nans::Bool`  which
determines if `NaNs` should be added in book keeping.
`eigmethod` which can be `:eig`, `:eigs` or
`:invpow`. The `:invpow` is an implementation of the power method, which
is slow but works well e.g. for `BigFloat`. The kwarg `recompute_U`
determines if the `U`-matrix should be recomputed in every
deflation (which can be more robust). The implementation
has two loggers `logger` and `inner_logger`. The `logger`
corresponds to outer iterations (deflation)
and  `inner_logger` is the iterates in Broydens method.
The kwarg `check_error_every` and `print_error_every`
  detemine how often errors should be check and how
often they should be printed. For real problems
with complex conjugate symmetry, you may want
to set the kwarg `addconj=true` in order
to reduce computation by automatically adding the
complex conjugate vectors.


The method computes an invariant pair and can therefore find
several eigenvalues. The
retured value is (S,V) is an invariant pair and
the eigenvalues are on the diagonal of S.
Eigenpairs can be directly extracted with
[`get_deflated_eigpairs`](@ref).



# Example

```julia-repl
julia> nep=nep_gallery("dep0");
julia> S,V=broyden(nep);
julia> λ=S[1,1]
-0.15955391823299253 - 3.874865266487398e-19im
julia> minimum(svdvals(compute_Mder(nep,λ)))
1.6293996560844023e-16
julia> λ=S[2,2]
-0.5032087003825461 + 1.1969823800738464im
julia> minimum(svdvals(compute_Mder(nep,λ)))
1.1073470346550144e-15
julia> λ=S[3,3]
1.2699713558173726 + 5.342786996459857e-16im
julia> minimum(svdvals(compute_Mder(nep,λ)))
5.905315846211231e-16
julia> broyden(nep,logger=2,check_error_every=1);  # Prints out a lot more convergence info
```
In order to extract eigenpairs you can use the following:
```julia-repl
julia> (D,X)=get_deflated_eigpairs(S,V,size(nep,1));
julia> for i=1:3; @show norm(compute_Mlincomb(nep,D[i],X[:,i])); end
norm(compute_Mlincomb(nep, D[i], X[:, i])) = 8.459878994614521e-13
norm(compute_Mlincomb(nep, D[i], X[:, i])) = 1.2102336671048442e-13
norm(compute_Mlincomb(nep, D[i], X[:, i])) = 2.1012363973403225e-16
```
# References

* Jarlebring, Broyden’s method for nonlinear eigenproblems, SIAM J. Sci. Comput., 41:A989–A1012, 2019, https://arxiv.org/pdf/1802.07322

"""
broyden(nep::NEP;params...)=broyden(ComplexF64,nep,:eye;params...)
broyden(nep::NEP,approxnep;params...)=broyden(ComplexF64,nep,approxnep;params...)
broyden(T,nep::NEP;params...)=broyden(T,nep,:eye;params...)
function broyden(::Type{TT},nep::NEP,approxnep::Union{NEP,Symbol,Matrix};σ::Number=0,
                 pmax::Integer=3,
                 c::Vector=ones(TT,size(nep,1)),
                 maxit::Integer=1000,addconj=false,
                 check_error_every::Integer=10,
                 print_error_every::Integer=1,
                 threshold::Real=0.2,
                 tol::Real=1e-12,
                 errmeasure::Function=broyden_default_errmeasure,
                 add_nans::Bool=false,
                 include_restart_timing::Bool=true,
                 eigmethod::Symbol=:eig,
                 logger =0,
                 recompute_U = false,
                 inner_logger = 0
                 ) where {TT<:Number}

    @parse_logger_param!(logger)
    @parse_logger_param!(inner_logger)

    time0=time_ns();
    n=size(nep,1);
    if (pmax>size(nep,1))
        @warn "Too many eigenvalues requested. Reducing"
        pmax=size(nep,1);
    end
    σ=TT(σ);


    # Step 1. Compute M0 and T0
    if (isa(approxnep,Matrix))
        M1=approxnep;
    elseif (approxnep == :eye)
        M1=Matrix{TT}(I,n,n)
    else
        M1=compute_Mder(approxnep,σ);
    end

    T1=inv(M1);


    X=zeros(TT,n,0);
    S=zeros(TT,0,0);

    k=1;

    all_errhist=[]; sumiter=1;
    all_timehist=[];
    all_iterhist=[];
    sumiter=1;
    UU = Matrix{TT}(I, n, pmax+1) # For storage
    U1=view(UU,1:n,1:0);
    while (k<=pmax)



        # Step 5
        ## Complete recomputation of U0
        p_U1=size(U1,2);
        U1=view(UU,1:n,1:k-1);
        if recompute_U
            i_start=1;
        else
            i_start=p_U1+1; # In theory they are already satisfied
        end
        for i=i_start:k-1
            ei=zeros(TT,size(S,1)); ei[i]=1;
            f = (σ * Matrix(1.0I, k-1, k-1) - S) \ ei
            vv=X*f;
            U1[:,i]=compute_Mlincomb(nep,σ,vv);
        end

        # Step 6
        MM::Matrix{TT}=[M1 U1; X' zeros(TT,k-1,k-1)];
        push_info!(logger,"running eigval comp for deflation",continues=true);

        local d,V;
        if (eigmethod==:eig)
            d,V = eigen(MM)
        elseif (eigmethod==:eigs)
            d,V = eigs(MM,which=:SM)
        elseif (eigmethod==:invpow)
            d,V = eigs_invpow(MM,maxit=4000,sigma=0)
        else
            error("Unknown eig method",eignmethod)
        end

        push_info!(logger,"");

        x=V[:,argmin(abs.(d))];


        if (TT <: Real)
            x=real(x);
        end


        # Orthogonalize
        v0=x[1:n];
        u0=x[n+1:end];
        h=X'*v0;
        v0=v0-X*h;
        u0 = u0 + (σ * Matrix{TT}(I, k-1, k-1) - S) * h
        CH=[X';c'];

        # Normalize it
        u0=u0/(c'*v0);
        v0=v0/(c'*v0);


        if (!include_restart_timing)
            time0=time_ns();
        end

        # Step 7
        d=sqrt(eps(real(TT)));
        push_info!(logger,"Computing initial matrix",continues=true);
        f1a=(compute_Mlincomb(nep,σ+d,v0)-compute_Mlincomb(nep,σ-d,v0))/2d;
        f1b = -U1 * ((σ * Matrix{TT}(I, k-1, k-1) - S) \ u0)
        push_info!(logger,".");
        f1=f1a+f1b;
        W1=[U1 f1];


        push_info!(logger,"Starting broyden n=$n")
        T=copy(T1);

        # Main call to broyden_T
        (λm,vm,um,Tm,Wm,iter,errhist,timehist)=
        broyden_T(TT,nep,
                  v1=v0,u1=u0,λ1=σ,
                  CH=CH,T1=T1,W1=W1,
                  S=S,X=X,
                  maxit=maxit,
                  check_error_every=check_error_every,
                  print_error_every=print_error_every,
                  threshold=threshold,
                  tol=tol,
                  errmeasure=errmeasure,
                  time0=time0,
                  logger=inner_logger)

        if size(all_iterhist,1) > 0
            iterhist=(1:size(errhist,1)) .+ all_iterhist[end]
        else
            iterhist=(1:size(errhist,1))
        end

        # Deflation book keeping

        if (add_nans && size(all_iterhist,1)>1)
            all_errhist=[all_errhist;NaN]
            all_timehist=[all_timehist;NaN]
            all_iterhist=[all_iterhist;NaN]
        end
        all_errhist=[all_errhist;errhist];
        all_timehist=[all_timehist;timehist];
        all_iterhist=[all_iterhist;iterhist];
        sumiter=sumiter+iter;
        um=um/norm(vm[1:n])
        vm=vm/norm(vm[1:n]) # Normalize

        push_info!(logger,"Found an eigval $k:$λm");


        X=[X vm]
        S=[S um;zeros(1,k-1) λm]


        # Add conjugate?
        if (abs(imag(λm))>tol*10 && addconj)


            v1 = conj(vm + X[:,1:k-1] * ((λm * Matrix{TT}(I, k-1, k-1) - S[1:k-1,1:k-1]) \ um))
            λ1=conj(λm);

            rnorm=norm(compute_Mlincomb(nep,λ1,v1))
            push_info!(logger,"Adding conjugate $k")
            if (rnorm>tol*10)
                @warn "Trying to add a conjugate pair which does not have a very small residual."
            end

            h=X'*v1;
            v1t=v1-X*h;
            beta=norm(v1t);
            X=[X v1t/beta];
            k=k+1;
            S1=zeros(TT,k,k);
            S1[1:(k-1),1:(k-1)]=S;
            S1[k,k]=λ1;
            R = Matrix{TT}(I, k, k)
            R[1:k-1,end]=h; R[k,k]=beta;
            S=(R*S1)/R;
            #println("opnorm(XX-I)=",opnorm(X'*X-Matrix(1.0I, size(X,2), size(X,2))))
            #println("opnorm(MM)=",opnorm(compute_MM(nep,S,X)));
            #        X=[X v1]
            #      S=[S conj(v1[(n0+1):(n0+p)]);zeros(1,p) conj(λ1)]
            #dnep=NEPBroydenDeflated(nep,S,X);
        end
        k=k+1;
    end
    push_info!(logger,"Iterations:$sumiter")
    return S,X,T1,all_errhist,all_timehist,all_iterhist;

end





function eigs_invpow(MM;maxit=10,sigma=0)
    AA = factorize(MM - sigma * sparse(1.0I, size(MM,1), size(MM,1)))
    z=ones(size(AA,1));
    for k=1:maxit
        z=AA\z;
        z=z/opnorm(z);
    end
    lambda=z'*MM*z;
    z=reshape(z,size(AA,1),1);
    return ([lambda],z);
end
