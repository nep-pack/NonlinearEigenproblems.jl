module Gallery_old
  using NEPCore
  using PolynomialRoots
  using SparseArrays
  export nep_gallery
  """
  Returns a NEP object from a gallery of examples of nonlinear eigenvalue problems. name decides which NEP. \\
  Usage:\\
     nep=nep_gallery(name)\\
     or
     nep=nep_gallery(name,params)
  """
  function nep_gallery(name,params...)
      if (name == "pep0")
          # A polynomial eigenvalue problem

          n=200; # mat size
          p=3; # Poly degree

          A=Array{Float64}(n,n,p)
          Random.seed!(0)
          for i=1:p
              A[:,:,i]=randn(n,n);
          end

          # Create the derivative function for the PEP
          function PEP_Md(λ,i=0)
              # Only workds for i=0 or i=1
              if (i==0)
                  M=zeros(n,n)
                  for i=1:p
                      M+=λ^(i-1)*A[:,:,i]
                  end
                  return M
              elseif (i==1)
                  Mp=zeros(n,n)
                  for i=2:p
                      Mp += (i-1)*(λ^(i-2)*A[:,:,i])
                  end
                  return Mp
              else
                  error("PEP higher derivatives not yet implemented")
              end

          end

          nep=NEP(n,PEP_Md);

          nep.rf=function (x; y=x,
                           target=0,
                           λ0=target)
              c=zeros(p)
              for k=1:p
                  c[k]=dot(y,A[:,:,k]*x);
              end
              r=roots(c);
              x,index=findmin(abs(r-target))
              return r[index]
          end
          return nep
      elseif (name == "dep0")
          # A delay eigenvalue problem
          n=5;
          Random.seed!(0) # reset the random seed
          A0=randn(n,n);
          A1=randn(n,n);
          I=eye(n,n);
          tau=1;

          # Derivative function for DEPs
          DEP_Md=function DEP_Md(λ,i=0)
              if (i==0)
                  return -λ*I+A0+A1*exp(-tau*λ)
              elseif (i==1)
                  return -I-tau*A1*exp(-tau*λ)
              else
                  return ((-tau)^i)*A1*exp(-tau*λ)
              end
          end
          nep=NEP(n,DEP_Md);

          return nep
      elseif (name == "dep0_sparse")
          # A delay eigenvalue problem with sparse matrices
          n=5;
          Random.seed!(0) # reset the random seed
          A0=sparse(randn(n,n));
          A1=sparse(randn(n,n));
          I=sparse(eye(n,n));
          tau=1;

          # Derivative function for DEPs
          function DEP_Md_sparse(λ,i=0)
              if (i==0)
                  return -λ*I+A0+A1*exp(-tau*λ)
              elseif (i==1)
                  return -I-tau*A1*exp(-tau*λ)
              else
                  return ((-tau)^i)*A1*exp(-tau*λ)
              end
          end
          nep=NEP(n,DEP_Md_sparse);

          return nep
      else
          error("NEP with name '"*name*"' not found in gallery.")
      end
  end
end
