module NEPSolver
  using NEPCore
  export newtonraphson
  
  function newtonraphson(nep::NEP;
                         errmeasure=NaN,
                         tolerance=1e-10,
                         maxit=10,
                         λ=0,
                         v=randn(nep.n,1),
                         c=v,
                         displaylevel=0)
      if (isnan(errmeasure))
          # If no relresnorm available use resnorm
          if (isdefined(nep, :relresnorm))
              errmeasure=nep.relresnorm;
          else
              errmeasure=nep.resnorm;
              if (displaylevel>0)
                  println("Using resnorm")
              end
          end
      end
      
      for k=1:maxit
          if (displaylevel>0)
              println("Iteration:",k,
                      " resnorm:",
                      errmeasure(λ,v))
          end
          if (errmeasure(λ,v)<tolerance)
              return (λ,v)
          end

          # Compute NEP matrix and derivative 
          M=nep.Md(λ)
          Md=nep.Md(λ,1)

          # Create jacobian
          J=[M Md*v; c' 0];
          F=[M*v; c'*v-1];

          # Compute update
          delta=-J\F;

          # Update eigenvalue and eigvec
          v=v+delta[1:nep.n];
          λ=λ+delta[nep.n+1];
              
      end

      # TODO: throw No convergence exception
      
  end

      
end

