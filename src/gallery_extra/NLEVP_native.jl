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
