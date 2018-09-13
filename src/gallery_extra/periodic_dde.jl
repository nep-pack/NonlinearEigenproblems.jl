using SparseArrays
using Random

import Base.size;
import ..NEPCore.compute_Mder
import ..NEPCore.compute_Mlincomb
import ..NEPCore.compute_Mlincomb!
import ..NEPCore.compute_MM
import ..NEPCore.compute_Mlincomb_from_MM
import ..NEPCore.compute_Mlincomb_from_Mder




"""
   type PeriodicDDE_NEP <: NEP

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
abstract type PeriodicDDE_NEP <: NEP; end;
struct PeriodicDDE_NEP_DAE <: PeriodicDDE_NEP
    A::Function
    B::Function
    n::Integer
    N::Integer
    tau::Number;
    disconts::Vector{Float64};
    E::Matrix;
    isconst::Bool
    function PeriodicDDE_NEP_DAE(A::Function,B::Function,E::AbstractMatrix,tau::Number)
        n=size(A(0),1)
        nep=new(A,B,n,1000,tau,[],E,false);;
        return nep;
    end
end


struct PeriodicDDE_NEP_ODE <: PeriodicDDE_NEP
    A::Function
    B::Function
    n::Integer
    N::Integer
    tau::Number;
    disconts::Vector{Float64};
    #ode_solver::Function
    E::Matrix
    function PeriodicDDE_NEP_ODE(A::Function,B::Function,tau::Number)
        n=size(A(0),1)
        nep=new(A,B,n,1000,tau,[]);;
        #nep.disconts=detect_disconts(nep);
        return nep;
    end
end



struct NEP_Mder <: NEP
    Mder::Function
    n::Number
end
function compute_Mder(nep::NEP_Mder,λ::Number,der::Integer=0)
    print("Computing Mder");
    MM=nep.Mder(λ);
    println(".");
    return MM
end

function size(nep::NEP_Mder)
    n = nep.n
    return (n,n)
end
function size(nep::NEP_Mder,dim)
    n = nep.n
    return n
end

function size(nep::PeriodicDDE_NEP)
    n = nep.n
    return (n,n)
end
function size(nep::PeriodicDDE_NEP,dim)
    n = nep.n
    return n
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
        y[:]=y+(s1+2*s2+2*s3+s4)/6;
        t=t+h;
    end
    return y;
end

function ode_be_dae(A,E,a,b,N,y0)
    h=(b-a)/N;
    y=copy(y0);
    t=a+h;
    for k=1:N
        y=(h*A(t)-E)\(E*y);
        t=t+h;
    end
    return y;
end

#function ode_be_dae_const(A,E,a,b,N,y0)
#    h=(b-a)/N;
#    y=copy(y0);
#    t=a+h;
#    AA=A(t);
#    Afact=factorize(h*AA-E);
#    for k=1:N
#        y=Afact\(E*y);
#        t=t+h;
#    end
#    return y;
#end
#
#abstract type ODESolver end;
#type RK4Solver <: ODESolver; end
#
#function forward_solve(nep::PeriodicDDE_NEP,S,V,Type{RK4Solver})
#
#    h=(b-a)/N; t=a;
#    y=copy(y0);
#    for k=1:N
#        s1=h*f(t,y);
#        s2=h*f(t+h/2,y+s1/2);
#        s3=h*f(t+h/2,y+s2/2);
#        s4=h*f(t+h,y+s3);
#        y=y+(s1+2*s2+2*s3+s4)/6;
#        t=t+h;
#    end
#    return y;
#
#end
#


#function compute_Mlincomb(nep::PeriodicDDE_NEP, λ::Number, v::Any)
#    n=size(nep,1);
#    # We are using (non-trivial) fact that
#    # the MM satisfies an ODE (as well as the action)
#    F=(t,Y) -> (nep.A(t)*Y+nep.B(t)*Y*exp(-(nep.tau*λ))-Y*λ)
#    Y0=v;
#    YY=ode_rk4(F, 0,nep.tau, nep.N, Y0);
#    #YY=myode(F, 0,nep.tau, nep.N, Y0);
#    return YY-Y0
#end
#

function compute_MM(nep::PeriodicDDE_NEP_DAE, S ,V)
    if (size(V,2)>1)
        println(size(V));
        error("Not implemented");
    end
    Af= t -> (nep.A(t)+nep.B(t)*exp(-nep.tau*S[1,1])-S[1,1]*nep.E);
    Y0=V;
    if (!nep.isconst)
        YY=ode_be_dae(Af,nep.E, 0,nep.tau, nep.N, Y0);
        return YY-Y0;
    else
        YY=ode_be_dae_const(Af,nep.E, 0,nep.tau, nep.N, Y0);
        return YY-Y0;
    end

end


### The implementations of compute_MM is a
### crude first-version, not optimized for efficiency nor accuracy.
function compute_MM(nep::PeriodicDDE_NEP_ODE, S ,V)
    n=size(nep,1);
    # We are using (non-trivial) fact that
    # the MM satisfies an ODE (as well as the action)
    local F::Function;
    if size(S,1)==1
        F=(t,Y) -> (nep.A(t)*Y+nep.B(t)*Y*exp(-nep.tau*S[1,1])-Y*S[1,1])
    else
        F=(t,Y) -> (nep.A(t)*Y+nep.B(t)*Y*exp(-Matrix(nep.tau*S))-Y*S)
    end

    Y0=V; # Use RK4 for ODE's
    YY=ode_rk4(F, 0,nep.tau, nep.N, Y0);
    #YY=myode(F, 0,nep.tau, nep.N, Y0);
    return YY-Y0
end

function compute_Mlincomb(nep::PeriodicDDE_NEP,λ::Number,V::Union{AbstractVector,AbstractMatrix})
    return compute_Mlincomb_from_MM(nep,λ,V,ones(eltype(V),size(V,2)))
end


##
#function compute_MM(nep::PeriodicDDE_NEP, S ,V)
#    n=size(nep,1);
#    # We are using (non-trivial) fact that
#    # the MM satisfies an ODE (as well as the action)
#    F=(t,Y) -> (nep.A(t)*Y+nep.B(t)*Y*exp(-full(nep.tau*S))-Y*S)
#    Y0=V;
#    t0=0;
#    disconts_and_tau=[nep.disconts;nep.tau];
#    YY=copy(V);
#    for k=1:size(disconts_and_tau,1)
#        YY=ode_rk4(F, t0,disconts_and_tau[k], nep.N, Y0);
#        Y0=YY;
#        t0=disconts_and_tau[k];
#    end
#    #YY=myode(F, 0,nep.tau, nep.N, Y0);
#    return YY-V
#end


# For compute_Mlincomb we (implicitly) use compute_Mlincomb_from_MM

function compute_Mder(nep::PeriodicDDE_NEP,λ::Number,der::Integer=0)
    if (der==0)
        n = size(nep,1)
        ek = zeros(n); ek[1] = 1
        z1 = compute_Mlincomb(nep,λ,ek)
        Z = zeros(eltype(z1), n, n)
        Z[:,1] = z1
        for k=2:n
            ek=copy(Z[:,k]); ek[k]=1
            Z[:,k]=compute_Mlincomb(nep,λ,ek)
        end
        return Z
        #return compute_Mder_from_MM(nep,λ,der)
    elseif (der==1)
        # Compute first derivative with finite difference. (Slow and inaccurate)
        ee=sqrt(eps())/10
        Yp=compute_Mder(nep,λ+ee,0)
        Ym=compute_Mder(nep,λ-ee,0)
        return (Yp-Ym)/(2ee)
    else # This would be too slow.
        error("Higher derivatives not implemented");
    end

end



"""
    nep=problem_gallery(PeriodicDDE_NEP; name)

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
function periodic_dde_gallery(::Type{PeriodicDDE_NEP}; name::String="mathieu",n=200,TT=ComplexF64)
    if (name == "mathieu")
        δ=1; b=1/2; a=0.1; tau=2;
        A=t-> [0 1; -( δ+ a*cos(pi*t) ) -1];
        B=t-> [0 0; b 0];
        nep=PeriodicDDE_NEP_ODE(A,B,tau)
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
        A0=sprandn(n,n,0.3)-one(TT)*I
        A1=sprandn(n,n,0.3)-one(TT)*I
        B0=sprandn(n,n,0.3)-one(TT)*I
        B1=sprandn(n,n,0.3)-one(TT)*I
        tau=2;
        A=t-> A0+cos(pi*t)*A1;
        B=t-> B0+exp(0.01*sin(pi*t))*B1;
        nep=PeriodicDDE_NEP_ODE(A,B,tau)
        return nep;
    elseif (name == "discont")
        δ=1; b=1/2; a=0.1; tau=2;
        A=t-> [0 1; -( δ+ a*cos(pi*t) ) -1]+Matrix(1.0*I,2,2)*((t-0.3)^2)*(t>0.3);
        B=t-> [0 0; b 0];
        nep=PeriodicDDE_NEP_ODE(A,B,tau)
    elseif (name == "milling0")
        δ=1; b=1/2; a=0.1; tau=2;
        A0=t-> [0 1; -( δ+ a*cos(pi*t) ) -1]+Matrix(1.0*I,2,2)*((t-0.3)^2)*(t>0.3);
        B0=t-> [0 0; b 0];
        e=ones(n-1);
        DD=spdiagm(-1 => e, 0 => -2*e, 1 => e)
        DD[1,1]=-1;
        h=1/n;
        DD=DD
        a=1;
        A=t-> [sparse(A0(t)) a*speye(2,n); a*speye(n,2) DD];
        B=t-> [sparse(B0(t)) spzeros(2,n); spzeros(n,2) -0*speye(n)];
        #DD=full(spdiagm(-1 => e, 0 => -2*e, 1 => e))
        #A=t-> [(A0(t)) eye(2,n); eye(n,2) DD];
        #B=t-> [(B0(t)) zeros(2,n); zeros(n,2) eye(n)];

        nep=PeriodicDDE_NEP(A,B,1)
    elseif (name == "milling1_be")
        # The milling model in analyzed by Insperger, Orosz, Bueler, etc (with unit constants)
        omega0=TT(1);
        zeta0=TT(1);  # Zeta (damping)
        m  =TT(1);      # Mass
        ap =TT(1)
        KR =TT(1)
        KT =TT(1);
        tau=TT(1)

        A0=[TT(0) TT(1); -omega0^2 -2*zeta0*omega0];

        phi=t -> 2*pi*t/tau;

        h=t -> (t<real(tau)/2).*(sin(phi(t)).^2*KR+KT*cos(phi(t)).*sin(phi(t)));
        E21=zeros(TT,2,2); E21[2,1]=1;

        # We need to use implicit method (default for DAEs).
        nep=PeriodicDDE_NEP_DAE(t-> A0-E21*h(t)*ap/m,
                                t->+E21*h(t)*ap/m,Matrix(one(TT)*I,2,2),1)

        nep.N=50;
        return nep;
    elseif (name == "milling1_rk4")
        # The milling model in analyzed by Insperger, Orosz, Bueler, etc (with unit constants)
        omega0=TT(1);
        zeta0=TT(1);  # Zeta (damping)
        m  =TT(1);      # Mass
        ap =TT(1)
        KR =TT(1)
        KT =TT(1);
        tau=TT(1)

        A0=[TT(0) TT(1); -omega0^2 -2*zeta0*omega0];

        phi=t -> 2*pi*t/tau;

        h=t -> (t<tau/2).*(sin(phi(t)).^2*KR+KT*cos(phi(t)).*sin(phi(t)));
        E21=zeros(TT,2,2); E21[2,1]=1;

        # We need to use implicit method (default for DAEs).
        nep=PeriodicDDE_NEP_ODE(t-> A0-E21*h(t)*ap/m,
                                t->+E21*h(t)*ap/m,1)
        nep.N=50;
        return nep;



    elseif (name == "milling")
        # The problem by Rott and Hömberg (and Jarlebring)
        #error("Problem in http://dx.doi.org/10.3182/20100607-3-CZ-4010.00023 not yet implemented")
        m=setDefaultParameters(n); m=computeMatricesScaled(m);
        ap=0.1;nn=3200;

        #x=1e-2;d=[ones(m.disc_n0+1);x*ones(m.disc_n0+1)];
        #D=eye(2*m.disc_n0+2);


        A=t-> -AScaled(t,nn/60,ap,m)/100;
        B=t-> -BScaled(t,nn/60,ap,m)/100;
        nep=PeriodicDDE_NEP_ODE(A,B,1)
        return nep;
    elseif (name == "newmilling")
        error("Not yet implemented");
    else
        error("Unknown PeriodicDDE_NEP type:",name);
    end
end
