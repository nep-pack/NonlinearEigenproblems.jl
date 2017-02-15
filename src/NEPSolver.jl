module NEPSolver
    using NEPCore
    using NEPTypes
    using MATLAB
    using LinSolvers
    ## NEP-Methods
    export compsolve #Wrapper around the solver for a linearized PEP pencil

    include("method_newton.jl")
    include("method_iar.jl")

    #############################################################################
    #Solve the linearized companion of a PEP
    function compsolve(pep::PEP)

        #Linearize to Ax = λEx
        E,A = companion(pep);

        #Choose eigensolver
        local eigsolverfunc::Function;
        
        levsolver = LinEigSolver();

        if (issparse(pep))
            eigsolverfunc = matlab_eigs_PEP;
        else
            eigsolverfunc = julia_eig_PEP;
        end

        #D,V = eigsolverfunc(A,E);
        D,V = levsolver.solve(A,B=E,λ_t=1.0);

        return D,V
    end

    #############################################################################
    #Call MATLAB eigs() for the linearized PEP
    function matlab_eigs_PEP(A,E)
        AA=mxarray(A)
        EE=mxarray(E)
        num = mxarray(size(A)[1]);

        @mput AA EE num
        @matlab begin
            n = num(1);

            (V,D) = eigs(AA,EE,n);
        end
        @mget D V

        return diag(D),V
    end

    #############################################################################
    #Call Julia eig() for the linearized PEP 
    function julia_eig_PEP(A,E)
        D,V = eig(A,E);
        return D,V;
    end

    ### Moved to NEPCore.jl
    ##############################################################################
    #  function default_errmeasure(nep::NEP, displaylevel)
    #      # If no relresnorm available use resnorm
    #      if (isdefined(nep, :relresnorm))
    #          return nep.relresnorm;
    #      else
    #          if (displaylevel>0)
    #              println("Using resnorm")
    #          end
    #          return nep.resnorm;
    #      end
    #  end
    #
        



end #End module
