module NEPSolver
    using NEPCore
    using NEPTypes
    using LinSolvers

    ## NEP-Methods
   
    include("method_newton.jl")
    include("method_iar.jl")
    include("method_mslp.jl")
    include("method_companion.jl")


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
