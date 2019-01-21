#  A delay eigenvalue problem with one eigenvalue equal to one.
function dep1()
    A0=([1 2 3 ; 4 5 6; 1 -1 3]);
    A1=((-A0+[1 0 3;0 0 -1;0 0 10])*exp(1));
    Q=[1 0 3; 1 1 -4; 2 3 1];
    A0=Q\(A0*Q);
    A1=Q\(A1*Q);
    nep=DEP([A0,A1],[0,1])
    return nep;
end


# A symmetric delay eigenvalue problem with double eigenvalues
# Examle from H. Voss and M. M. Betcke, Restarting iterative projection methods for Hermitian nonlinear eigenvalue problems with minmax property, Numer. Math., 2017
function dep_symm_double(n::Integer=100)
    LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

    x = range(0, stop = pi, length = n)
    h=x[2]-x[1];
    h=pi
    LL=LL/(h^2)
    LL=kron(LL,LL)

    b=broadcast((x,y)->100*abs(sin(x+y)),x,transpose(x))
    a=broadcast((x,y)->-8*sin(x)*sin(y),x,transpose(x))
    B=sparse(1:n^2,1:n^2,b[:])
    A=LL+sparse(1:n^2,1:n^2,a[:])

    nep=DEP([A,B],[0,2.0])
    return nep
end


# A delay eigenvalue problem with a double non-semisimple eigenvalue in λ=3πi
# Examle from E. Jarlebring, Convergence factors of Newton methods for nonlinear eigenvalue problems, LAA, 2012
function dep_double()
    denom = 8+5*pi;
    a1 = 2/5 *(65*pi + 32)/(denom);
    a2 = 9*pi^2*(13+5*pi)/(denom);
    a3 = 324/5 *pi^2*(5*pi+4)/(denom);
    b1 = (260*pi + 128 + 225*pi^2)/(10*denom);
    b2 = 45*pi^2/denom;
    b3 = 81*pi^2*(40*pi + 32 + 25*pi^2)/(10*denom);
    A0 = [ 0    1    0;  0    0    1;  -a3  -a2  -a1];
    A1 = [ 0    0    0;  0    0    0;  -b3  -b2  -b1];

    tau=1.0;
    nep=DEP([A0,A1],[0,tau])
    return nep
end


#  A quadratic problem with real eigenvalues
# Four smallest eigenvalues of the problem:
# 1.0e+03 * [ -2.051741417993845  -0.182101627437811  -0.039344930222838  -0.004039879577113]
function real_quadratic()
    A0 = [ 4     0     1     1;
           0     2     1     1;
           1     1     6    -2;
           1     1    -2     3];
    A1 = [ 167  -140    95  -131;
          -140   327    54    85;
            95    54   235   -81;
          -131    85   -81   181];
    A2 = [ 2     1    -1    -1;
           1     5    -3     2;
          -1    -3     3     0;
          -1     2     0     3];
    A=[A0,A1,A2]
    nep=PEP(A)
    return nep
end


# A quadratic delay eigenvalue problem in "The infinite Bi-Lanczos method for nonlinear eigenvalue problems",  Sarah W. Gaaf and Elias Jarlebring
function qdep0()
    qdepbase=joinpath(dirname(@__FILE__()),
                      "converted_misc" ,"qdep_infbilanczos_")
    A0=read_sparse_matrix(qdepbase * "A0.txt")
    A1=read_sparse_matrix(qdepbase * "A1.txt")
    tau::Float64 = 1
    quadfun = S -> S^2
    constfun = S -> one(S)
    expfun = S -> exp(-tau*S)

    AA = [-one(A0), A0, A1]
    fi = [quadfun, constfun, expfun]
    return SPMF_NEP(AA, fi)
end


# A quadratic delay eigenvalue problem in "A linear eigenvalue algorithm for the  nonlinear eigenvalue problem",  Elias Jarlebring, Wim Michiels, Karl Meerbergen \\
function qdep1()
    A0=[ 0.3000   -0.6000         0    0.4000
        -0.3000    0.4000   -0.8000    1.9000
         0.1000   -1.6000   -1.3000         0
        -1.4000   -0.9000    0.2000    0.9000];
    A1=[ 0.8000    0.2000   -1.3000   -0.3000
        -1.1000    0.9000    1.2000    0.5000
         0.5000    0.2000   -1.6000   -1.3000
         0.7000    0.4000   -0.4000         0];
    return SPMF_NEP([one(A0), A0, A1], [λ -> -λ^2, λ -> one(λ), λ -> exp(-λ)])
end


# A quadratic eigenvalue problem with chosen eigenvalues
function qep_fixed_eig(n::Integer=5, E::AbstractVecOrMat=NaN*ones(2*n))
    if any(isnan.(E))
        Random.seed!(0) # reset the random seed
        E=randn(2*n)
    end
    A1 = diagm(0 => E[1:n])
    A2 = diagm(0 => E[n+1:2*n])
    K = one(A1)
    nep=PEP([A1*A2,-A1-A2,K])
    return nep
end


# This problem stems from
# L. P. Shayer and S. A. Campbell.  Stability, bifurcation and multistability in a system of two coupled neurons with multiple time delays. SIAM J. Applied Mathematics , 61(2):673–700, 2000
# It is also a benchmark example in DDE-BIFTOOL
function neuron0()
    pars = [1/2; -1; 1; 2.34; 0.2; 0.2 ; 1.5] .+ 0im
    kappa = pars[1]
    beta = pars[2]
    A = [0 pars[3]; pars[4] 0]

    x = [0; 0]  # The zero (trivial) stationary solution
    # A non-trivial stationary solution
    #x=[3.201081590416643561697725111745656884148241428177442574927999582405266342752249e-01
    #   5.096324796647208606096018689631125587762848405395086474417800152349531876959548e-01]

    tauv = [0;0.2;0.2;1.5]
    A0 = -kappa * Matrix(1.0I, 2, 2)
    A1 = A[2,1] * [0 0; (1-tanh(x[2])^2) 0]
    A2 = A[1,2] * [0 (1-tanh(x[1])^2); 0 0]
    A3 = beta * diagm(0 => [(1-tanh(x[1])^2), (1-tanh(x[2])^2)])
    return DEP([A0, A1, A2, A3], tauv)
end


# DEP modelling a beam
function beam(n::Integer=100)
    h=1/n;
    ee = ones(n);
    A0 = spdiagm(-1 => ee[1:n-1], 0 => -2*ee, 1 => ee[1:n-1]);
    A0[end,end]=1/h;
    A0[end,end-1]=-1/h;
    A1=sparse([n],[n],[1.0]); # A1=en*en'
    tau=1.0;
    return DEP([A0,A1],[0,tau]);
end


function sine_nep()
    data_dir=joinpath(dirname(@__FILE__()),  "converted_sine")
    A0=read_sparse_matrix(joinpath(data_dir,"sine_A0.txt"));
    A1=read_sparse_matrix(joinpath(data_dir,"sine_A1.txt"));
    A2=read_sparse_matrix(joinpath(data_dir,"sine_A2.txt"));
    V=Matrix(read_sparse_matrix(joinpath(data_dir,"sine_V.txt")));
    Q=Matrix(read_sparse_matrix(joinpath(data_dir,"sine_Q.txt")));

    n=size(A0,1);
    Z=spzeros(n,n);
    pep=PEP([A0,A1,Z,Z,A2]);
    # Matrix  sine function. Note that the Term is rank two which is not exploited here
    sin_nep=SPMF_NEP([V*Q'], [S-> sin(S)]);
    nep=SPMFSumNEP(pep,sin_nep) # Note: nep has a low-rank term
    return nep;
end


# From the tutorial on NEP-PACK online user's manual (except we use LowRankNEP)
function schrodinger_movebc(n=1000,L0=1,L1=8,α=25*pi/2,V0=10.0)

    # Discreatization matrices
    xv=Vector(range(0,stop=L0,length=n))
    h=xv[2]-xv[1];
    n=size(xv,1);
    V=x->1+sin(α*x);
    Dn=spdiagm(-1 => [ones(n-2);0]/h^2, 0 => -2*ones(n-1)/h^2, 1 => ones(n-1)/h^2)
    Vn=spdiagm(0 => [V.(xv[1:end-1]);0]);
    In=spdiagm(0 => [ones(n-1);0])
    # Construct functions
    f1=S->one(S);
    f2=S->-S;
    hh=S-> sqrt(S+V0*one(S))
    g=S-> cosh((L1-L0)*hh(S))
    f=S-> inv(hh(S))*sinh((L1-L0)*hh(S))

    # First NEP
    nep1=SPMF_NEP([Dn-Vn,In],[f1,f2]);

    # Create the Low-Rank NEP
    Lv1=sparse([zeros(n-1,1);1.0]);
    Lv2=sparse([zeros(n-1,1);1.0]);
    Uv1=sparse([zeros(n-1,1);1.0]);
    Uv2=sparse([zeros(n-3,1); [1/(2*h), -2/h, 3/(2*h)]])
    nep2=LowRankFactorizedNEP([Lv1,Lv2],[Uv1,Uv2],[g,f]);

    # Create a sum of the two
    nep=SumNEP(nep1,nep2);
    return nep
end
