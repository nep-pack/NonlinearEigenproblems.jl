# Methods / functions in common for several contour integral methods
using Distributed, LinearAlgebra, Random


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
