using LinearAlgebra
#using Arpack # disabled. Use ArnoldMethod instead
using ..NEPCore, ..NEPTypes
import ..NEPCore.compute_Mder;
import ..NEPCore.compute_Mlincomb;
import Base.size;
export broyden;

abstract type NEPBroydenDeflated <: NEP  end

struct NEPBroydenDeflatedEll1 <: NEPBroydenDeflated;
    orgnep::NEP;
    S::AbstractArray;
    X::AbstractArray;
    function NEPBroydenDeflatedEll1(orgnep,S::AbstractArray,X)
        return new(orgnep,reshape(S,size(S,1),size(S,1)),
                   reshape(X,size(orgnep,1),size(X,2)))
    end
end

struct NEPBroydenDeflatedEll2 <: NEPBroydenDeflated;
    orgnep::NEP;
    S::AbstractArray;
    X::AbstractArray;
    function NEPBroydenDeflatedEll2(orgnep,S::AbstractArray,X)
        return new(orgnep,reshape(S,size(S,1),size(S,1)),
                   reshape(X,size(orgnep,1),size(X,2)))
    end
end

function compute_Mder(nep::NEPBroydenDeflated,λ::Number,i::Integer=0)
    if (i>0)
        error("Not implemented");
    end
    M=zeros(typeof(λ),size(nep,1),size(nep,1));
    for k=1:size(nep,1);
        ek=zeros(typeof(λ),size(nep,1));
        ek[k]=1;
        M[:,k]=compute_Mlincomb(nep,λ,ek);
    end
    return M
end

function deflated_errmeasure(nep::NEP,λ,v)
    return norm(compute_Mlincomb(nep,λ,v))/norm(v);
end

function extract_eigenpair(nep::NEPBroydenDeflatedEll1,λ)
    S=nep.S
    X=nep.X
    dd,VV = eigen(S)
    I=argmin(abs.(dd-λ))
    w=X*VV[:,I]; w=w/norm(w);
    return (λ,w);
end

function deflated_errmeasure(nep::NEPBroydenDeflatedEll1,λ,v)
    n0=size(nep.orgnep,1);
    p=size(nep.X,2);
    if (p==0)
        return norm(compute_Mlincomb(nep.orgnep,λ,v))/norm(v)
    end

    S=[nep.S v[(n0+1):(n0+p)];zeros(1,p) λ]
    X=[nep.X v[1:n0]];
    dd,VV = eigen(S)
    I=argmin(abs.(dd-λ))
    w=X*VV[:,I];

    return norm(compute_Mlincomb(nep.orgnep,λ,w))/norm(w)
end

function size(nep::NEPBroydenDeflated,dim=-1)
    n0=size(nep.orgnep,1);
    n=n0+size(nep.S,1);
    if (dim==-1)
        return (n,n)
    else
        return n
    end
end

function compute_Mlincomb(nep::NEPBroydenDeflatedEll1,λ::Number,x::AbstractArray)
    n0=size(nep.orgnep,1);
    p=size(nep.S,1);
    b1=x[1:n0];
    b2=x[(n0+1):(n0+p)];
    local z;
    if (p==0)
        z=b1
        f=compute_Mlincomb(nep.orgnep,λ,z);
    else
        z = b1 + nep.X*((λ*Matrix{eltype(λ)}(I, p, p) - nep.S) \ b2)
        f=vcat(compute_Mlincomb(nep.orgnep,λ,z),nep.X'*b1);
    end
    return Array{eltype(x),1}(f);
end

function compute_Mlincomb(nep::NEPBroydenDeflatedEll2,λ::Number,x::AbstractArray)
    n0=size(nep.orgnep,1);
    p=size(nep.S,1);
    b1=x[1:n0];
    b2=x[(n0+1):(n0+p)];
    local z;
    if (p==0)
        z=b1
        f=compute_Mlincomb(nep.orgnep,λ,z);
    else
        z = b1 + nep.X*((λ*Matrix{eltype(λ)}(I, p, p) - nep.S) \ b2)
        f1=compute_Mlincomb(nep.orgnep,λ,z);
        #f2=nep.X'*b1; # ell1
        # ell2
        f2=nep.X'*b1+λ*nep.S'*nep.X'*b1; # ell1
        f2=f2+nep.S'*(nep.X'*(nep.X*b2));
        f=vcat(f1,f2);
    end
    return Array{eltype(x),1}(f);
end

function broyden_naive_H(::Type{TT},nep::NEPBroydenDeflated;
                         v1=0,u1=[],λ1=TT(0),
                         CH=0,T1=0,W1=0,
                         maxit=100,
                         S=zeros(TT,0,0),X=zeros(TT,0,0),
                         check_error_every=10,
                         print_error_every=1,
                         tol=1e-12,
                         threshold=0.4,
                         time0=time_ns(),
                         errmeasure::Function=broyden_default_errmeasure,
                         logger=0
                         ) where {TT<:Number}

    @parse_logger_param!(logger)

    n=size(nep.orgnep,1);
    p=size(nep,1)-n;

    T1=Array{TT,2}(T1);
    CH=Array{TT,2}(CH);
    v1=Array{TT,1}(v1);
    λ1=TT(λ1);

    M1=inv(T1);
    c=zeros(TT,n+p);
    c[1:n]=CH[end:end,:]';
    #K=[M1 W1[1:n,1:p]; CH[1:p,:] zeros(TT,p,p)];
    #k=zeros(TT,n+p);
    #k[1:n]=W1[:,p+1];
    J=[M1 W1;CH zeros(TT,p+1,p+1)]
    H=inv(J);


    x=Array{TT,1}([v1;u1;λ1]);
    F=vcat(compute_Mlincomb(nep,x[n+p+1],x[1:(n+p)]),0)
    #(y1,y2)=fake_Mlincomb(nep,x[n+p+1],x[1:n],x[n+1:n+p]);
    #F=vcat(y1,y2);
    errhist=Array{real(TT),1}(NaN*ones(real(TT),maxit));
    timehist=Vector{Float64}(maxit)*NaN;
    II = Matrix{TT}(I, p, p)

    for j=1:maxit

        Δx=-H*F; # called s_k


        γ=TT(1.0);
        tt=norm(Δx);
        if (tt>threshold)  # Avoid too big λ-steps
            γ=TT(threshold)/tt;
        end
        #        γ=TT(1.0);
        #        if (abs(Δx[end])>threshold)  # Avoid too big λ-steps
        #            γ=TT(threshold)/abs(Δx[end]);
        #        end

        xp=x+γ*Δx;

        # Impose structure => This improves performance (at least for deflation)
        #if (impose_normalization)
        #    xp[1:n]=xp[1:(n+p)]/(c'*xp[1:(n+p)]); Δx=xp-x;
        #end

        #(y1,y2)=fake_Mlincomb(nep,xp[n+p+1],xp[1:n],xp[n+1:n+p]);
        #Fp=vcat(y1,y2);
        Fp=vcat(compute_Mlincomb(nep,xp[n+p+1],xp[1:(n+p)]),c'*xp[1:(n+p)]-1)
        ΔF=Fp-F;

        # Update Jacobian approximation
        z1=(Fp-(1-γ)*F)/γ;
        Hz1=H*z1;
        aH=(Δx'*H)/(norm(Δx)^2+Δx'*Hz1);
        H=H-Hz1*aH

        # step forward
        x=xp;
        F=Fp;

        if (mod(j,check_error_every)==0)
            #errhist[j]=deflated_errmeasure(nep,x[n+p+1],x[1:n+p])
            #errhist[j]=opnorm(Fp);
            λ=x[end];  vv=x[1:nep.orgnep.n];  uu=x[(nep.orgnep.n+1):end-1];
            errhist[j]=errmeasure(λ,vv+nep.X*((λ*II-nep.S)\uu),  F);
            timehist[j]=Float64((time_ns()-time0)*1e-9);
            if (mod(j,print_error_every)==0)
                d = opnorm(CH*x[1:n] - reverse(Matrix{TT}(I, 1+p, 1), dims = 1))
            end

            #println(j," Normf=",opnorm(F), " λ=",xp[n+p+1]);
            if (errhist[j]<tol)
                return (x[n+p+1],x[1:n],x[(n+1):(n+p)],H,0,j,errhist[1:j],timehist[1:j])
            end

        end

    end

    push_info!(logger,"Too many iterations"); #" resnorm=",opnorm(rk));
    #error("Too many iterations")
    return (x[n+p+1],x[1:n],x[(n+1):(n+p)],H,H[1:n,(n+1):end],maxit,errhist[1:end])


end


function broyden_naive_J(::Type{TT},nep::NEPBroydenDeflated;
                         v1=0,u1=[],λ1=TT(0),
                         CH=0,T1=0,W1=0,
                         maxit=100,
                         S=zeros(TT,0,0),X=zeros(TT,0,0),
                         check_error_every=10,
                         print_error_every=1,
                         tol=1e-12,
                         threshold=0.4,
                         time0=time_ns(),
                         errmeasure::Function=broyden_default_errmeasure,
                         logger=0
                         ) where {TT<:Number}

    @parse_logger_param!(logger)

    n=size(nep.orgnep,1);
    p=size(nep,1)-n;

    T1=Array{TT,2}(T1);
    CH=Array{TT,2}(CH);
    v1=Array{TT,1}(v1);
    λ1=TT(λ1);
    u1=Array{TT,1}(u1);

    M1=inv(T1);
    c=zeros(TT,n+p);
    c[1:n]=CH[end:end,:]';
    #K=[M1 W1[1:n,1:p]; CH[1:p,:] zeros(TT,p,p)];
    #k=zeros(TT,n+p);
    #k[1:n]=W1[:,p+1];
    J=[M1 W1;CH zeros(TT,p+1,p+1)];
    if (isa(nep,NEPBroydenDeflatedEll2))
        J[n+1:end-1,:]=[nep.X'+λ1*nep.S'*nep.X'  nep.S'*nep.X'*nep.X nep.S'*nep.X'*v1];
    end


    x=Array{TT,1}([v1;u1;λ1]);
    F=vcat(compute_Mlincomb(nep,x[n+p+1],x[1:(n+p)]),TT(0))
    errhist=Array{real(TT),1}(NaN*ones(real(TT),maxit));
    timehist=Vector{Float64}(maxit)*NaN;
    II = Matrix{TT}(I, p, p)



    for j=1:maxit

        Δx=-J\F; # called s_k


        γ=TT(1.0);
        tt=norm(Δx);
        if (tt>threshold)  # Avoid too big λ-steps
            γ=TT(threshold)/tt;
        end

        #γ=TT(1.0);
        #if (abs(Δx[end])>threshold)  # Avoid too big λ-steps
        #    γ=TT(threshold)/abs(Δx[end]);
        #end

        xp=x+γ*Δx;


        # Impose structure => This improves performance (at least for deflation)
        #if (impose_normalization)
        #    xp[1:n]=xp[1:(n+p)]/(c'*xp[1:(n+p)]); Δx=xp-x;
        #end

        Fp=vcat(compute_Mlincomb(nep,xp[n+p+1],xp[1:(n+p)]),c'*xp[1:(n+p)]-1)
        ΔF=Fp-F;

        # Update Jacobian approximation
        β=1/(norm(Δx)^2);


        #z1=(K*Δv+k*Δλ -ΔF[1:n]); # Can be simplified
        z1=(Fp-(1-γ)*F)/γ;


        J=J+z1*(Δx'*β);


        # step forward
        x=xp;
        F=Fp;

        if (mod(j,check_error_every)==0)
            λ=x[end];  vv=x[1:nep.orgnep.n];  uu=x[(nep.orgnep.n+1):end-1];
            errhist[j]=errmeasure(λ,vv+nep.X*((λ*II-nep.S)\uu),  F);
            #errhist[j]=deflated_errmeasure(nep,x[n+p+1],x[1:n+p])
            timehist[j]=Float64((time_ns()-time0)*1e-9);
            if (mod(j,print_error_every)==0)
                d = opnorm(CH*x[1:n] - reverse(Matrix{TT}(I, 1+p, 1), dims = 1))
            end
            #println(j," Normf=",opnorm(F), " λ=",xp[n+p+1]);
            if (errhist[j]<tol)
                return (x[n+p+1],x[1:n],x[(n+1):(n+p)],J,0,j,errhist[1:j],timehist[1:j])
            end

        end

    end

    push_info!(logger,"Too many iterations"); #" resnorm=",opnorm(rk));
    #error("Too many iterations")
    return (x[n+p+1],x[1:n],x[(n+1):(n+p)],J[1:n,1:n],J[1:n,(n+1):end],maxit,errhist[1:end])


end

function broyden_default_errmeasure(λ,v,r)
    return norm(r)/norm(v)
end


function broyden_ell2(::Type{TT},nep::NEP;
                      v1=0,u1=[],λ1=TT(0),
                      CH=0,T1=0,W1=0,
                      maxit=100,
                      S=zeros(TT,0,0),X=zeros(TT,size(nep,1),0),
                      check_error_every=10,
                      print_error_every=1,
                      tol=1e-12,
                      threshold=0.4,
                      time0=time_ns(),
                      errmeasure::Function=broyden_default_errmeasure) where {TT<:Number}
    # Check types are consistent

    v=Array{TT,1}(v1);
    u=Array{TT,1}(u1);
    X=Array{TT,2}(X);
    S=Array{TT,2}(S);
    CH=Array{TT,2}(CH);
    λ=TT(λ1);


end


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

        #        if (abs(Δλ[end])>threshold)  # Avoid too big λ-steps
        #            γ=TT(threshold)/abs(Δλ);
        #        end
        v0=v; u0=u; λ0=λ;

        v=v+γ*Δv;
        u=u+γ*Δu;
        λ=λ+γ*Δλ;


        # Impose structure => This improves performance (at least for deflation)
        #if (impose_normalization)
        #    h=X'*v;
        #    v=v-X*h;
        #    α=CH[end,:].'*v[1:n];
        #    v=v/α;
        #    u=u/α;
        #end

        # Update Jacobian approximation
        vv=v+(X*((λ*II-S)\u));
        rkp=compute_Mlincomb(nep,λ,vv);


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
            #errhist[j]=opnorm(rkp);
            #errhist[j]=opnorm(rk)/opnorm(v+(X/(λ*II-S))*u);
            if (mod(j,print_error_every)==0)
                d = opnorm(CH*v - reverse(Matrix{TT}(I, 1+p, 1), dims = 1))
            end

            #if (d>1e-10)
            #    println(CH*v - reverse(Matrix{TT}(I, 1+p, 1), dims = 1))
            #end


            push_iteration_info!(logger,j,err=errhist[j]);
            if (errhist[j]<tol)
                return (λ,v,u,T,W,j,errhist[1:j],timehist[1:j])
            end
        end

    end

    push_info!(logger,"Too many iterations"); #" resnorm=",opnorm(rk));
    #error("Too many iterations")
    #return (λ,[v;u],zeros(n+1,n+1),j,errhist[1:j])
    return (λ,v,u,T,W,maxit,errhist,timehist)



end

##################
"""
    S,V = broyden([eltype,]nep::NEP[,approxnep::NEP];kwargs)

Runs Broydens method (with deflation) for the nonlinear eigenvalue problem defined by nep.
An approximate nep can be provided which is used as an initialization of starting
matrix/vectors.

The method computes an invariant pair and can therefore find several eigenvalues. The
retured value is (S,V) is an invariant pair and the eigenvalues are on the diagonal of S.

See [`newton`](@ref) for other parameters.


# Example

```julia-repl
julia> nep=nep_gallery("dep0");
julia> S,V=broyden(nep);
julia> λ=S[1,1]
-0.3587189459686267 - 3.0010731412746105e-31im
julia> minimum(svdvals(compute_Mder(nep,λ)))
1.6066157878930856e-16
julia> λ=S[2,2]
-0.04093521177097334 + 1.486011530941621im
julia> minimum(svdvals(compute_Mder(nep,λ)))
4.159109513753696e-16
julia> λ=S[3,3]
0.8347353572199486 + 1.5032076225139986e-14im
julia> minimum(svdvals(compute_Mder(nep,λ)))
1.296144276122994e-14
julia> broyden(nep,logger=2,check_error_every=1);  % Prints out a lot more convergence info
```

# References

* Jarlebring, Broyden’s method for nonlinear eigenproblems, 2018, https://arxiv.org/pdf/1802.07322

"""
broyden(nep::NEP;params...)=broyden(ComplexF64,nep,nep;params...)
broyden(nep::NEP,approxnep::NEP;params...)=broyden(ComplexF64,nep,approxnep;params...)
function broyden(::Type{TT},nep::NEP,approxnep::NEP;σ::Number=0,
                 pmax::Integer=3,
                 c::Vector=ones(TT,size(nep,1)),
                 maxit::Integer=1000,addconj=false,
                 check_error_every::Integer=10,
                 print_error_every::Integer=1,
                 broyden_variant::Symbol=:T,
                 threshold::Real=0.2,
                 tol::Real=1e-12,
                 errmeasure::Function=broyden_default_errmeasure,
                 add_nans::Bool=false,
                 include_restart_timing::Bool=true,
                 eigmethod::Symbol=:eig,
                 logger =0,
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

    M1=compute_Mder(approxnep,σ);
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
        for i=(p_U1+1):k-1
            #for i=1:k-1 # If you want to recompute
            ei=zeros(TT,size(S,1)); ei[i]=1;
            f = (σ * Matrix(1.0I, k-1, k-1) - S) \ ei
            vv=X*f;
            U1[:,i]=compute_Mlincomb(nep,σ,vv);
        end

        # Step 6
        MM::Matrix{TT}=[M1 U1; X' zeros(TT,k-1,k-1)];
        push_info!(logger,"running eig",continues=true);

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


        # Not in MS yet:
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


        if (broyden_variant == :T)
            push_info!(logger,"Running T variant *********************************** n=$n")
            T=copy(T1);

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
        elseif (broyden_variant == :J)
            push_info!(logger,"Running J variant *********************************** n=$n")
            dnep=NEPBroydenDeflatedEll1(nep,S,X);

            (λm,vm,um,Tm,Wm,iter,errhist,timehist)=
            broyden_naive_J(TT,dnep,
                            v1=v0,u1=u0,λ1=σ,
                            CH=CH,T1=T1,W1=W1,
                            S=S,X=X,
                            maxit=maxit,
                            check_error_every=check_error_every,
                            print_error_every=print_error_every,
                            threshold=threshold,
                            tol=tol,
                            time0=time0,
                            logger=inner_logger)
        elseif (broyden_variant == :H)
            push_info!(logger,"Running H variant *********************************** n=$n")
            dnep=NEPBroydenDeflatedEll1(nep,S,X);

            (λm,vm,um,Tm,Wm,iter,errhist,timehist)=
            broyden_naive_H(TT,dnep,
                            v1=v0,u1=u0,λ1=σ,
                            CH=CH,T1=T1,W1=W1,
                            S=S,X=X,
                            maxit=maxit,
                            check_error_every=check_error_every,
                            print_error_every=print_error_every,
                            threshold=threshold,
                            tol=tol,
                            time0=time0,
                            logger=inner_logger)
        else
            error("Unknown broyden method");
        end

        if size(all_iterhist,1) > 0
            iterhist=(1:size(errhist,1)) .+ all_iterhist[end]
        else
            iterhist=(1:size(errhist,1))
        end
        add_nans=true;


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
        #println("Quality of eigval guess:", abs(λ0-λ1)/abs(λ1))
        #I=argmin(abs.(λv-λ1))
        #println("Best guess distance:", abs(λv[I]-λ1)/abs(λ1))
        #println("It was ",λv[I])
        #println(" not   ",λ0)


        X=[X vm]
        S=[S um;zeros(1,k-1) λm]

        #println("J=",[inv(Tm) Wm; CH zeros(k,k)]);

        #println("opnorm(MM)=",opnorm(compute_MM(nep,S,X)));
#println("I-X'*X=",opnorm(Matrix(1.0I, k, k)-X'*X))

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



function deflated_broyden_ell2(::Type{TT},nep::NEP,approxnep::NEP;σ=0,
                               pmax::Integer=3,
                               c=ones(TT,size(nep,1)),
                               maxit=1000,addconj=false,
                               check_error_every=10,
                               print_error_every=1,
                               broyden_variant=:J,threshold=0.2,
                               tol=1e-12,
                               errmeasure::Function=broyden_default_errmeasure,
                               add_nans=false,
                               include_restart_timing=true
                               ) where {TT<:Number}

    time0=time_ns();
    n=size(nep,1);
    if (pmax>2*size(nep,1))
        @warn "Too many eigenvalues requested. Reducing"
        pmax=size(nep,1);
    end
    σ=ComplexF64(σ);


    # Step 1. Compute M0 and T0

    M1=compute_Mder(approxnep,σ);
    T1=inv(M1);


    X=zeros(ComplexF64,n,0);
    S=zeros(ComplexF64,0,0);

    k=1;

    all_errhist=[]; sumiter=1;
    all_timehist=[];
    all_iterhist=[];
    sumiter=1;
    UU = Matrix{TT}(I, n, pmax+1) # For storage
    U1=view(UU,1:n,1:0);
    while (k<=pmax)



        # Step 5
        p_U1=size(U1,2); # Not documented
        U1=view(UU,1:n,1:k-1);
        for i=(p_U1+1):k-1
            #for i=1:k-1
            ei=zeros(TT,size(S,1)); ei[i]=1;
            U1[:,i] = compute_Mlincomb(nep, σ, X * ((σ*Matrix{TT}(I, k-1, k-1) - S) \ ei))
        end

        # Step 6
        MM=[M1 U1; X'+σ*S'*X' S'*(X'*X)];
        d,V = eigen(MM)
        x=V[:,argmin(abs.(d))];

        # Not in MS yet:
        # Orthogonalize
        v0=x[1:n];
        u0=x[n+1:end];
        #h=X'*v0;
        #v0=v0-X*h
        CH=[X';c'];

        # Normalize it
        u0=u0/(c'*v0);
        v0=v0/(c'*v0);


        if (!include_restart_timing)
            time0=time_ns();

        end

        # Step 7
        d=sqrt(eps(real(TT)));
        f1a=(compute_Mlincomb(approxnep,σ+d,v0)-compute_Mlincomb(approxnep,σ-d,v0))/2d;
        f1b = -U1 * ((σ * Matrix{TT}(I, k-1, k-1) - S) \ u0)
        f1=f1a+f1b;
        W1=[U1 f1];


        if (broyden_variant == :J)
            println("Running J variant *********************************** n=",n);
            dnep=NEPBroydenDeflatedEll2(nep,S,X);


            (λm,vm,um,Tm,Wm,iter,errhist,timehist)=broyden_naive_J(ComplexF64,dnep,
                                                                   v1=v0,u1=u0,λ1=σ,
                                                                   CH=CH,T1=T1,W1=W1,
                                                                   S=S,X=X,
                                                                   maxit=maxit,
                                                                   check_error_every=check_error_every,
                                                                   print_error_every=print_error_every,                                              threshold=threshold,
                                                                   tol=tol,
                                                                   time0=time0,
                                                                   logger=inner_logger)

        elseif (broyden_variant == :H)
            println("Running H variant *********************************** n=",n);
            dnep=NEPBroydenDeflatedEll2(nep,S,X);


            (λm,vm,um,Tm,Wm,iter,errhist,timehist)=broyden_naive_H(ComplexF64,dnep,
                                                                   v1=v0,u1=u0,λ1=σ,
                                                                   CH=CH,T1=T1,W1=W1,
                                                                   S=S,X=X,
                                                                   maxit=maxit,
                                                                   check_error_every=check_error_every,
                                                                   print_error_every=print_error_every,                                              threshold=threshold,
                                                                   tol=tol,
                                                                   time0=time0,
                                                                   logger=inner_logger)

        else
            error("Unknown broyden method");
        end

        #println("Should_be_zero=",[X X*S]'*[X
        if (size(all_iterhist,1)>0)
            iterhist=(1:size(errhist,1))+all_iterhist[end];
        else
            iterhist=(1:size(errhist,1))
        end
        add_nans=true;


        if (add_nans && size(all_iterhist,1)>1)
            all_errhist=[all_errhist;NaN]
            all_timehist=[all_timehist;NaN]
            all_iterhist=[all_iterhist;NaN]
        end
        all_errhist=[all_errhist;errhist];
        all_timehist=[all_timehist;timehist];
        all_iterhist=[all_iterhist;iterhist];
        sumiter=sumiter+iter;
        β=opnorm(vcat(vm,X*um+λm*vm));
        um=um/β
        vm=vm/β # Normalize
        println("Found an eigval ",k,":",λm);
        #println("Quality of eigval guess:", abs(λ0-λ1)/abs(λ1))
        #I=argmin(abs.(λv-λ1))
        #println("Best guess distance:", abs(λv[I]-λ1)/abs(λ1))
        #println("It was ",λv[I])
        #println(" not   ",λ0)


        X=[X vm]
        S=[S um;zeros(1,k-1) λm]

        println("normalize:",[X;X*S]'*[X;X*S])

        println("S=",S)

        #println("J=",[inv(Tm) Wm; CH zeros(k,k)]);

        #println("opnorm(MM)=",opnorm(compute_MM(nep,S,X)));
        #println("I-X'*X=",opnorm(Matrix(1.0I, k, k)-X'*X))

        if (abs(imag(λm))>sqrt(eps()) && addconj)


            v1 = conj(vm + X[:,1:k-1] * ((λm * Matrix{TT}(I, k-1, k-1) - S[1:k-1,1:k-1]) \ um))
            λ1=conj(λm);

            rnorm=norm(compute_Mlincomb(nep,λ1,v1))
            println("Adding conjugate ",k,
                    " norm(res)=",rnorm);
            if (rnorm>tol*10)
                @warn "Trying to add a conjugate pair which does not have a very small residual."
            end

            h=X'*v1;
            v1t=v1-X*h;
            beta=opnorm(v1t);
            X=[X v1t/beta];
            k=k+1;
            S1=zeros(ComplexF64,k,k);
            S1[1:(k-1),1:(k-1)]=S;
            S1[k,k]=λ1;
            R = Matrix{ComplexF64}(I, k, k)
            R[1:k-1,end]=h; R[k,k]=beta;
            S=(R*S1)/R;
            #println("opnorm(XX-I)=",opnorm(X'*X-Matrix(1.0I, size(X, 2), size(X, 2))))
            #println("opnorm(MM)=",opnorm(compute_MM(nep,S,X)));
            #        X=[X v1]
            #      S=[S conj(v1[(n0+1):(n0+p)]);zeros(1,p) conj(λ1)]
            #dnep=NEPBroydenDeflated(nep,S,X);
        end
k=k+1;
end
println("Iterations:",sumiter)
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
