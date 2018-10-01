export gauss_legendre_weights;   # For debugging only
export distributed_kernel_trapezoidal;    # For debugging only
export distributed_kernel_gauss_legendre; # For debugging only

"
#Creates the NEP associated with example in E. Jarlebring and
#W. Michiels and K. Meerbergen, The infinite  {Arnoldi} method and an application to time-delay systems with distributed delays, Delay Systems - Methods, Applications and New Trends, 2012
"
function gallery_dep_distributed()
    # E. Jarlebring and W. Michiels and K. Meerbergen},
    # The infinite  {Arnoldi} method and an application to time-delay systems with distributed delays,
    # Delay Systems - Methods, Applications and New Trends
    # 2012
    #
    # Some correct eigenvalues:
    # -0.400236388049641 + 0.970633098237807i
    # -0.400236388049641 - 0.970633098237807i
    #  2.726146249832675 + 0.000000000000000i
    # -1.955643591177653 + 3.364550574688863i
    # -1.955643591177653 - 3.364550574688863i
    #  4.493937056300693 + 0.000000000000000i
    # -1.631513006819252 + 4.555484848248613i
    # -1.631513006819252 - 4.555484848248613i
    # -1.677320660400946 + 7.496870451838560i
    # -1.677320660400946 - 7.496870451838560i
    #
    #



    A0 = Matrix(-1.0I, 3, 3)
    A1=[2.5    2.8   -0.5
        1.8    0.3    0.3
        -2.3   -1.4    3.5];
    A2=[1.7    0.7   -0.3
        -2.4   -2.1   -0.2
        2.0    0.7    0.4];
    A3=[1.4   -1.3    0.4
        1.4    0.7    1.0
        0.6    1.6    1.7];
    idop= S -> S
    oneop = S -> one(S)
    f1 = S -> exp(-S)
    N = 10
    f2 = S -> distributed_kernel_gauss_legendre(S, N)
    #N = 1000
    #f2 = S -> distributed_kernel_trapezoidal(Matrix(S), N)
    return SPMF_NEP([A0,A1,A2,A3],[idop,oneop,f1,f2])
end


function distributed_kernel_trapezoidal(S,N0)

    fS = x -> exp(x*S) * (exp((x+0.5)^2) - exp(1/4))

    F=zero(S)

    h=1/N0;
    for i=1:(N0+1) # Trapezoidal rule
        x=float((i-1)*h-1)
        if (i==1 || i==N0+1)
            F += 0.5*fS(x)*h;
        else
            F += fS(x)*h;
        end
    end
    return F
end
"""
Computes distributed kernel matrix function using
Gauss-Legender quadrature.
"""
function distributed_kernel_gauss_legendre(S,N)
    f=x-> (exp((x+0.5)^2)-exp(1/4))
    F=zero(S);
    xv,wv=gauss_legendre_weights(N,-1,0);
    local E = one(S)

    accumulative_expm_comp=true
    for i=1:length(xv)
        if (accumulative_expm_comp)
            # An accumulative way to compute E=exp(xv[i]*S) which is
            # faster due to the fact that scaling and squaring
            # for exp((xv[i]-xv[i-1])*S) is faster than
            # exp(xv[i]*S)
            if (i==1)
                E=exp(xv[1]*S);
            else
                E=E*exp((xv[i]-xv[i-1])*S)
            end
        else
            E=exp(xv[i]*S);

        end

        fSw=E*(f(xv[i])*wv[i]);
        F=F+fSw;
    end
    return F
end




function  gauss_legendre_weights(N,a,b)

    N1=N; N2=N+1;
    xu = range(-1, stop = 1, length = N1)

    # L will be the Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);



    # "Derivative" of L
    Lp=zeros(N1,N2);


    # Starting values of Newton's method
    y = cos.((2 * (0:(N-1)) .+ 1) * pi / (2 * (N-1) + 2)) + (0.27/N1) * sin.(pi * xu * (N-1) / N2)
    y0 = 2


    local Lp0
    # Newton's method to decide points
    # Iterate until new points are uniformly within epsilon of old points

    while maximum(abs.(y .- y0)) > eps()



        L[:,1]=ones(size(L,1));
        Lp[:,1]=zeros(size(Lp,1));

        L[:,2]=y;
        Lp[:,2]=ones(size(Lp,1));

        # Recurrence
        for k=2:N1
            L[:,k+1]=( (2*k-1)*(y).*L[:,k]-(k-1)*L[:,k-1] )/k;
        end

        Lp0 = N2 * (L[:,N1] - y .* L[:,N2]) ./ (1 .- y.^2)

        y0=y;
        y=y0-L[:,N2]./Lp0;
    end

    # Linear map from[-1,1] to [a,b]
    x = (a * (1 .- y) + b * (1 .+ y)) / 2

    # Compute the weights
    w = (b-a) ./ ((1 .- y.^2) .* Lp0.^2) * (N2/N1)^2

    return x,w
end
