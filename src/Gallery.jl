module Gallery
  using NEPCore
  using PolynomialRoots
  export nep_gallery
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
          n=5;
          srand(0) # reset the random seed
          A0=randn(n,n);
          A1=randn(n,n);
          I=eye(n,n);
          tau=1;
          nep=DEP([A0,A1],[0,tau])
          return nep
      elseif (name == "dep0_sparse")
          # A delay eigenvalue problem with sparse matrices
          n=5;
          srand(0) # reset the random seed
          A0=sparse(randn(n,n));
          A1=sparse(randn(n,n));
          I=sparse(eye(n,n));
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


      end
      
  end 
end
