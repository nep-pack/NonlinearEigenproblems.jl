# Time-periodic delay-differential equations
module GalleryPeriodicDDE
    using NonlinearEigenproblems.NEPCore
    using NonlinearEigenproblems.Gallery
    using SparseArrays
    using Random

    import Base.size
    import NonlinearEigenproblems.NEPCore.compute_Mlincomb
    import NonlinearEigenproblems.NEPCore.compute_Mder
    import NonlinearEigenproblems.NEPCore.compute_MM
    import NonlinearEigenproblems.Gallery.nep_gallery

    export PeriodicDDE_NEP
    export nep_gallery

    """
   struct PeriodicDDE_NEP <: NEP

This type represents NEP associated with the time-periodic delay-differential equation
```math
̇{x}(t)=A(t)x(t)+B(t)x(t-τ)
```
where `A(t)` and `B(t)` are periodic functions with period `τ`.

# Example
```julia-repl
julia> nep=PeriodictDDE_NEP((t) -> 3, (t) -> 4,1)
julia> λ,v=newton(nep,λ=3)
julia> compute_Mlincomb(nep,λ,v)
1-element Array{Complex{Float64},1}:
 0.0+0.0im
```
"""
    struct PeriodicDDE_NEP <: NEP
        A::Function
        B::Function
        n::Integer
        N::Integer
        tau::Number
        function PeriodicDDE_NEP(A::Function,B::Function,tau::Number)
            n=size(A(0),1)
            return new(A,B,n,1000,tau);
        end
    end


    function size(nep::PeriodicDDE_NEP,dim=-1)
        n=nep.n;
        if (dim==-1)
            return (n,n)
        else
            return n
        end
    end


    function myode(f,a,b,N,y0)
        h=(b-a)/N
        t=a;
        yy=copy(y0);
        for k=1:N
            yy=yy+h*f(t,yy)
            t=t+h;
        end
        return yy;
    end

    # Stand-alone implementation of RK4
    function ode_rk4(f,a,b,N,y0)
        h=(b-a)/N; t=a;
        y=copy(y0);
        for k=1:N
            s1=h*f(t,y);
            s2=h*f(t+h/2,y+s1/2);
            s3=h*f(t+h/2,y+s2/2);
            s4=h*f(t+h,y+s3);
            y=y+(s1+2*s2+2*s3+s4)/6;
            t=t+h;
        end
        return y;
    end

    # The implementations of compute_MM is a
    # crude first-version, not optimized for efficiency nor accuracy.
    function compute_MM(nep::PeriodicDDE_NEP, S ,V)
        n=size(nep,1);
        # We are using (non-trivial) fact that
        # the MM satisfies an ODE (as well as the action)
        F=(t,Y) -> (nep.A(t)*Y+nep.B(t)*Y*exp(-Matrix(nep.tau*S))-Y*S)
        Y0=V;
        YY=ode_rk4(F, 0,nep.tau,nep.N,Y0);
        return YY-Y0
    end

    # For compute_Mlincomb we use compute_Mlincomb_from_MM
    compute_Mlincomb(nep::PeriodicDDE_NEP,λ::Number,
                 V::Union{AbstractMatrix,AbstractVector},a::Vector=ones(eltype(V),size(V,2)))=
             compute_Mlincomb_from_MM(nep,λ,V,a)


    function compute_Mder(nep::PeriodicDDE_NEP,λ::Number,der::Integer=0)
        if (der==0)
            return compute_Mder_from_MM(nep,λ,der)
        elseif (der==1)
            # Compute first derivative with finite difference. (Slow and inaccurate)
            ee=sqrt(eps())/10;
            Yp=compute_Mder(nep,λ+ee,0)
            Ym=compute_Mder(nep,λ-ee,0)
            return (Yp-Ym)/(2ee);
        else # This would be too slow.
            error("Higher derivatives not implemented");
        end

    end


"""
    nep=nep_gallery(PeriodicDDE_NEP; name)

Constructs a PeriodicDDE object from a gallery. The default is "mathieu", which is
the delayed Mathieu equation in section 1 in reference.

Subset of eigenvalues of the default example:
 -0.24470143590830754
 -0.561610418452567 - 1.511169478595549im
 -0.561610418452567 + 1.511169478595549im
 -1.859617846506182 - 1.261010754174415im
 -1.859617846506182 + 1.261010754174415im
 -2.400415351524992 - 1.323058902477801im
 -2.400415351524992 + 1.323058902477801im


# Example
```julia-repl
julia> nep=nep_gallery(PeriodicDDE_NEP,name="mathieu")
julia> (λ,v)=newton(Float64,nep,λ=-0.2,v=[1;1])
(-0.24470143590830754, [-1.25527, 0.313456])
julia> exp(nep.tau*λ)  # Reported in Figure 2 with multipliers in reference
0.6129923199095035
```
# Reference
* E. Bueler, Error Bounds for Approximate Eigenvalues of Periodic-Coefficient Linear Delay Differential Equations, SIAM J. Numer. Anal., 45(6), 2510–2536

"""
    function nep_gallery(::Type{PeriodicDDE_NEP}; name::String="mathieu",n=200)
        if (name == "mathieu")
            δ=1; b=1/2; a=0.1; tau=2;
            A=t-> [0 1; -( δ+ a*cos(pi*t) ) -1];
            B=t-> [0 0; b 0];
            nep=PeriodicDDE_NEP(A,B,tau)
            return nep;
        elseif (name == "rand0")
            # Some eigenvalues for n=200:
            #  4.63633+1.10239im
            # 4.63633-1.10239im
            # 5.58214+4.03225im
            # 5.58214-4.03225im
            # 5.73989+0.732386im
            # 5.73989-0.732386im

            Random.seed!(0);
            A0=sprandn(n,n,0.3)-speye(n,n)
            A1=sprandn(n,n,0.3)-speye(n,n)
            B0=sprandn(n,n,0.3)-speye(n,n);
            B1=sprandn(n,n,0.3)-speye(n,n);
            tau=2;
            A=t-> A0+cos(pi*t)*A1;
            B=t-> B0+exp(0.01*sin(pi*t))*B1;
            nep=PeriodicDDE_NEP(A,B,tau)
            return nep;
        elseif (name == "milling")
            # The problem by Rott and Hömberg (and Jarlebring)
            error("Problem in http://dx.doi.org/10.3182/20100607-3-CZ-4010.00023 not yet implemented")
        else
            error("Unknown PeriodicDDE_NEP type");
        end
    end
end
