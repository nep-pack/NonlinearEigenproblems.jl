module NEPSolver
    using NEPCore
    using NEPTypes
    using LinSolvers

    ## NEP-Methods
   
    include("method_newton.jl")
    include("method_iar.jl")
    include("method_mslp.jl")
    include("method_companion.jl")
    include("method_nlar.jl")

"""
     Computes an eigenvector approximation from an
     eigenvalue approximation (with very little
     computational effort).
"""
    function compute_eigvec_from_eigval(nep::NEP,λ)
        compute_Mder(nep,λ)
        # Not sure how to do this
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
