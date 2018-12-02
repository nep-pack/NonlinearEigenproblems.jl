# This file contains code for native implementations of some of the
# benchmarks in NLEIGS

function nlevp_native_gun()
    gunbase=joinpath(dirname(@__FILE__()),
      "converted_nlevp", "gun_")
    K=read_sparse_matrix(gunbase * "K.txt")
    M=read_sparse_matrix(gunbase * "M.txt")
    W1=read_sparse_matrix(gunbase * "W1.txt")
    W2=read_sparse_matrix(gunbase * "W2.txt")
    # The gun problem is a sum of a PEP and a problem containing square roots.
    pep=PEP([K,-M])
    sqrt1op= S -> 1im*sqrt(S)
    sqrt2op= S -> 1im*sqrt(S-108.8774^2*one(S))
    sqrtnep=SPMF_NEP([W1,W2],[sqrt1op,sqrt2op])
    nep=SumNEP(pep,sqrtnep)
    return nep
end


function nlevp_native_cd_player()
    cdbase=joinpath(dirname(@__FILE__()),
      "converted_nlevp", "cd_player_")
    K = Matrix(read_sparse_matrix(cdbase * "K.txt"))
    C = Matrix(read_sparse_matrix(cdbase * "C.txt"))
    M = one(K)
    nep = PEP([K, C, M])
    return nep
end


function nlevp_native_fiber()
    # Since the bessel functions are not available as matrix functions
    # we rely on interpolation (of "denominator" and "numerator" separately)
     # Construct the complicated third function
    L::Float64=2400;
    besselkp= (m,z)->  - besselk(m-1,z) - m*besselk(m,z)/z;  # Derivative
    numer= x-> ((L+0.5)/L^2)*x/(besselk(1,ComplexF64(x))^2)
    denom= x-> 1/(besselkp(1,ComplexF64(x))*besselk(1,ComplexF64(x)));
    f_sqrt= x ->  numer(x)/denom(x)
    s3=λ -> f_sqrt(sqrt(λ)*L)  # This is the complicated function in "fiber"

    m=10;
    TT=Complex{BigFloat}; # Do interpolation very accurately
    ### Different interpolation approaches ####
    ## Disc
    #phi=2im*pi*range(0.0;length=m+1,stop=1);
    #interp_points=0.1*exp.(phi[1:end-1]); #
    ## Uniform on an interval
    interp_points=0.01.+range(0;length=m,stop=3);; # This choice is done with eigenvalue information in mind (x approx 2.027892679678583)
    interp_points=Vector{TT}(interp_points) # High precision

    # Do the interpolation on numerator and denominator
    (Newton_Matrix,fnum)=construct_newton_matrix(TT,numer,interp_points)
    (Newton_Matrix,fdenom)=construct_newton_matrix(TT,denom,interp_points)
    num_coeffs=Newton_Matrix\fnum;
    denom_coeffs=Newton_Matrix\fdenom;

    # Now go back to ComplexF64
    interp_points_c64=Vector{ComplexF64}(interp_points);
    num_coeffs_c64=Vector{ComplexF64}(num_coeffs);
    denom_coeffs_c64=Vector{ComplexF64}(denom_coeffs);

    # Interpolated numerator and denominators
    numer_new_c64=x->newton_eval(num_coeffs_c64,x,interp_points_c64)
    denom_new_c64=x->newton_eval(denom_coeffs_c64,x,interp_points_c64)

    f_new_c64=x->  denom_new_c64(x)\numer_new_c64(x)
    s3_new= λ -> f_new_c64(sqrt(λ)*L); # This is the new function. Works for matrices and functions
    n=2400;

    A2=sparse([n],[n],[1.0])  # This is a rank-one matrix!
    A1=one(A2)

    # Create the A0-matrix takes a bit more work
    eta_cl = 1.4969;  # "physical params"
    alpha=25; ell=1.1;
    gam=0.003; delta=0.01;

    k_cl = 2*pi*eta_cl/ell;
    n_c = 400; n = 6*n_c;
    r = (1:n+1)*delta;
    mm = 1;

    inc = (1:n_c);
    i_n = (n_c+1:n-1);
    e = ones(n_c);

    # Helper functions. Note C(r) is r-independent, according to NLEVP.
    C = sqrt.( (1 .- 2*gam*(inc/n_c).^alpha) / (1 - 2*gam) ) .- 1;
    eta0 = r-> (eta_cl .+ 1.4201*C)
    k = r -> 2*pi*eta0(r)/ell;

    # setup the vectors in the diagonal
    y1 = -2*e - mm^2*(e ./ inc.^2) + delta^2*(k(r[1:n_c]).^2 .- k_cl^2);
    e = ones(size(i_n));
    y2 = -2*e - mm^2*(e ./ i_n.^2);
    y = [y1; y2; -1 + 1/(2*n) - mm^2/n^2];
    i = 1:n-1;
    z = (i.+0.5) ./ sqrt.( i.*(i.+1) );

    A0=spdiagm(0  => y[1:n], 1=> z[1:n-1], -1 => z[1:n-1])

    # The functions
    f1=S-> one(S);
    f2=S-> -S;
    f3=s3_new; # Interpolated function
    # f3=s3 # Use this if you want the exact formula

    nep=SPMF_NEP([A0,A1,A2],[f1,f2,f3])

    return nep;
end


# Helper functions for Newton interpolation (mainly for nlevp_native_fiber)
function construct_newton_matrix(T,ff,interp_points)
    m=size(interp_points,1);
    Newton_Matrix=zeros(T,m,m);
    Newton_Matrix[:,1] .= 1
    for col=2:m
        for row=col:m
            Newton_Matrix[row,col]=Newton_Matrix[row,col-1]*(interp_points[row]-interp_points[col-1])
        end
    end
    f=zeros(T,m);
    for k=1:m
        f[k]=ff(interp_points[k]);
    end
    return Newton_Matrix,f
end
function newton_eval(coeffs,S,interp_points)  # This works for λ::Number and λ::Matrix
    F=coeffs[1]*one(S);
    prod=one(S);
    for k=2:size(coeffs,1)
        prod = prod*(S-interp_points[k-1]*one(S))
        F += prod*coeffs[k];
    end
    return F
end



function nlevp_native_hadeler(α=100,n=8)
    i=1:n
    I2 = ones(n)*i'  # matrix with constant columns 1,2,3,4...n
    # Identity matrix
    II = Matrix{typeof(α)}(I,n,n);
    # Matrices
    A0= α*II
    A2 = n*II+1 ./(I2 + I2');
    B = ((n+1) .- max.(I2',I2)) .* (i*i');
    # Functions
    fv= [S->-one(S) ; S->S^2 ; S->(exp(S)-one(S))]
    nep=SPMF_NEP([A0, A2, B], fv);
    return nep
end

function nlevp_native_pdde_stability(n=15)

    # default values
    a0 = 2;
    b0 = .3;
    a1 = -2;
    b1 = .2;
    a2 = -2;
    b2 = -.3;
    t1 = -pi/2;

    # Construct a DEP
    h = pi/(n+1);  x = (1:n)*h;
    e = ones(n);
    A0 = spdiagm(-1 => e[1:end-1], 0 => -2*e, 1 => e[1:end-1])/(h^2)
    A0 = A0 + spdiagm(0 => (a0 .+ b0*sin.(x)))
    A1 = spdiagm(0 => (a1 .+ b1*x.*(1 .- exp.(x .- pi))))
    A2 = spdiagm(0 => (a2 .+ b2*x.*(pi .- x)))

    # The stabiltiy of the DEP given by -λI+A0+A1exp(-t1 λ)+A2exp(- TT λ) can be
    # as a function of TT can be analyzed with a QEP
    II=complex(sparse(1.0*I,n,n))
    E = kron(II,A2);
    gamma = exp(1im*t1); gamma = gamma/abs(gamma);
    F = kron(II,A0-gamma*A1) + kron(A0+gamma*A1,II);

    # Compute a permutation vector (faster than multiplication)
    p = vec(copy(reshape(1:n^2,n,n)'));

    Av = [E[p,p],F,E];

    return PEP(Av)
end




function toeplitz(v::AbstractVector)
    # This can be replaced with the types in ToeplitzMatrices
    # once the package is more stable (e.g. matrix matrix products works in general).
    n=length(v);
    T=zeros(eltype(v),n,n);
    for i=1:n
        for j=1:(n-i+1)
            T[i,j+i-1]=v[j]
            T[j+i-1,i]=v[j]
        end
    end
    if (issparse(v))
        return SparseMatrixCSC(T)
    else
        return T;
    end

end

function nlevp_native_loaded_string(n=20,kappa=1,m=1)
    # See Gallery for documentation
    A0 = toeplitz(SparseVector([2.0*n;-n;zeros(n-2)]));
    A1 = zeros(n,n); A1[n,n] = n-A0[n,n];
    B0 = toeplitz(SparseVector([4/(6*n);1/(6*n);zeros(n-2)]));
    B1 = zeros(n,n); B1[n,n] = 2/(6*n)-B0[n,n]
    C = zeros(n,n); C[n,n] = kappa;
    σ=kappa/m;
    fv=[S->one(S), S->-S, S -> (S-σ*one(S))\S]

    spmf1=SPMF_NEP([A0,B0],fv[1:2]);
    spmf2=SPMF_NEP([A1,B1,C],fv);  # This can be a Low-Rank SPMF

    nep=SPMFSumNEP(spmf1,spmf2)
    return nep;

end
