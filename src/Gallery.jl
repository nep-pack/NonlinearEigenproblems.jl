module Gallery
    using NEPCore
    using NEPTypes
    using MATLAB
    using PolynomialRoots
    
    export nep_gallery



    include("gallery_extra/distributed_example.jl")
    include("gallery_extra/nlevp_interface.jl")
    include("gallery_extra/waveguide.jl")
    
  """
  Returns a NEP object from a gallery of examples of nonlinear eigenvalue problems. name decides which NEP. \\
  Usage:\\
     nep=nep_gallery(name)\\
     or
     nep=nep_gallery(name,params)
  """
  function nep_gallery(name,params...)
      if (name == "dep0")
          # A delay eigenvalue problem
          local n
          if (length(params)>0) 
              n=params[1]
          else
              n=5; # Default size if 
          end

          srand(0) # reset the random seed
          A0=randn(n,n);
          A1=randn(n,n);
          tau=1;
          nep=DEP([A0,A1],[0,tau])
          return nep

      elseif (name == "dep0_sparse")
          # A delay eigenvalue problem with sparse matrices
          n=5;
          srand(0) # reset the random seed
          A0=sparse(randn(n,n));
          A1=sparse(randn(n,n));
          tau=1;
          nep=DEP([A0,A1],[0,tau])
          return nep

      elseif (name == "dep_double")
          # A delay eigenvalue problem with a double non-semisimple eigenvalue in λ=3πi
          # Examle from E. Jarlebring, Convergence factors of Newton methods for nonlinear eigenvalue problems, LAA, 2012
          n=3;

          denom = 8+5*pi;
          a1 = 2/5 *(65*pi + 32)/(denom);
          a2 = 9*pi^2*(13+5*pi)/(denom);
          a3 = 324/5 *pi^2*(5*pi+4)/(denom);
          b1 = (260*pi + 128 + 225*pi^2)/(10*denom);
          b2 = 45*pi^2/denom;
          b3 = 81*pi^2*(40*pi + 32 + 25*pi^2)/(10*denom);
          A0 = [ 0    1    0;  0    0    1;  -a3  -a2  -a1];
          A1 = [ 0    0    0;  0    0    0;  -b3  -b2  -b1];

          tau=1;
          nep=DEP([A0,A1],[0,tau])
          return nep

      elseif (name== "pep0")
          # A polynomial eigenvalue problem
          n=200; # mat size

          srand(0)
          A0=randn(n,n)
          A1=randn(n,n)
          A2=randn(n,n)
          A=[A0,A1,A2]
          nep=PEP(A)
          return nep

       elseif (name== "pep0_sparse_003")
          # A sparse polynomial eigenvalue problem
          n=200; # mat size

          srand(0)
          A0=sprandn(n,n,0.03)
          A1=sprandn(n,n,0.03)
          A2=sprandn(n,n, 0.03)
          A=[A0,A1,A2]
          nep=PEP(A)
          return nep

      elseif (name== "real_quadratic")
          # Create a quadratic problem with real eigenvalues          
          n=4; # mat size
	   
	  A0 = [4     0     1     1;
    		0     2     1     1;
   		1     1     6    -2;
    		1     1    -2     3];


	A1 = [167  -140    95  -131;
	     -140   327    54    85;
 	       95    54   235   -81;
 	     -131    85   -81   181];


	A2 =  [2     1    -1    -1;
  	       1     5    -3     2;
  	      -1    -3     3     0;
	      -1     2     0     3];

	# Four smallest eigenvalues of the problem:
	# 1.0e+03 *
	# -2.051741417993845
  	# -0.182101627437811
  	# -0.039344930222838
  	# -0.004039879577113


          A=[A0,A1,A2]
          nep=PEP(A)
          return nep

      elseif (name== "dep_distributed")
          return gallery_dep_distributed();

     elseif (name== "waveguide")
          return gallery_waveguide(params...);
      end
      
  end
        

    
end
