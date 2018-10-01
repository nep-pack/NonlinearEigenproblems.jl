module NonlinearEigenproblems

using Reexport

"Include the specified module, then `@reexport` it."
macro include_reexport(file_name, module_name)
    :( include($file_name); @reexport using .$module_name )
end

# Re-export core functionality, to remove the need for "using" statements. Non-core
# functionality not re-exported below can be accessed by "using" those modules,
# for example "using NonlinearEigenproblems.RKHelper".
# Note: Order matters below, due to inter-module dependencies.
include(joinpath("utils", "Serialization.jl"))
@include_reexport "NEPCore.jl" NEPCore
@include_reexport "NEPTypes.jl" NEPTypes
@include_reexport "LinSolvers.jl" LinSolvers
include(joinpath("rk_helper", "RKHelper.jl"))
@include_reexport "NEPSolver.jl" NEPSolver
@include_reexport "Gallery.jl" Gallery

end
