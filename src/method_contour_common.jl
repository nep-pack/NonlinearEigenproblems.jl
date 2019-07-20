# Methods / functions in common for several contour integral methods
using Distributed, LinearAlgebra, Random


export MatrixIntegrator, MatrixTrapezoidal;

abstract type MatrixIntegrator ; end

abstract type MatrixTrapezoidal <: MatrixIntegrator; end



function integrate_interval(ST::Type{MatrixTrapezoidal},::Type{T},f,gv,a,b,N,logger) where {T<:Number}
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    f1=f(t[1]);
    m=size(gv,1);
    S = zeros(T,size(f1)...,m)

    push_info!(logger,string(ST)*": computing G",continues=true)
    # Matrix of gv evaluations
    G = zeros(T,N,m);
    for i=1:m
        push_info!(logger,".",continues=true);
        push_info!(logger,2,"t=$t",continues=true);
        gg=gv[i].(Vector(t));
        G[:,i] = gg
    end
    push_info!(logger,".",continues=false);

    # Do the sum
    push_info!(logger,string(ST)*": summing terms",continues=true);
    for i = 1:N
        push_info!(logger,".",continues=true);
        if (i==1) # Already computed
            temp = f1;
        else
            temp = f(t[i])
        end
        for j=1:m
            S[:,:,j] += temp*G[i,j];
        end
    end
    push_info!(logger,"")
    return S * h
end



# The following integration methods have input and output:

# Returns (S0,S1) where S0 and S1 are certain integrals
# S0 approx int_a^b f(t) dt
# S1 approx int_a^b f(t)*g(t) dt
# Typically f(t) would return a vector or a matrix


#  Carries out Gauss quadrature (with N) discretization points
#  by call to @parallel
function quadg_parallel(f,g,a,b,N)
    x,w=gauss(N);
    # Rescale
    w=w*(b-a)/2;
    t=a+((x+1)/2)*(b-a);
    # Sum it all together f(t[1])*w[1]+f(t[2])*w[2]...
    S = @distributed (+) for i = 1:N
        temp = f(t[i])*w[i]
        [temp,temp*g(t[i])]
    end
    return S[1], S[2];
end

function quadg(f,g,a,b,N)
    x,w=gauss(N);
    # Rescale
    w=w*(b-a)/2;
    t=a+((x+1)/2)*(b-a);
    S0 = zero(f(t[1])); S1 = zero(S0)
    # Sum it all together f(t[1])*w[1]+f(t[2])*w[2]...
    for i = 1:N
        temp = f(t[i])*w[i]
        S0 += temp
        S1 += temp*g(t[i])
    end
    return S0, S1;
end


# Trapezoidal rule for a periodic function f
function ptrapz(f,g,a,b,N)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S0 = zero(f(t[1])); S1 = zero(S0)
    for i = 1:N
        temp = f(t[i])
        S0 += temp
        S1 += temp*g(t[i])
    end
    return h*S0, h*S1;
end


function ptrapz_parallel(f,g,a,b,N)
    h = (b-a)/N
    t = range(a, stop = b-h, length = N)
    S = @distributed (+) for i = 1:N
        temp = f(t[i])
        [temp,temp*g(t[i])]
    end
    return h*S[1], h*S[2];
end
