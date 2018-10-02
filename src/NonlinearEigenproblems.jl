module NonlinearEigenproblems

using Reexport

"Include the specified module, then `@reexport` it."
macro include_export(file_name, module_name)
    :( include($file_name); @reexport using .$module_name )
end

# Re-export core functionality, to remove the need for "using" statements.
# Non-core functionality not re-exported below can be accessed by "using"
# those modules, for example "using NonlinearEigenproblems.RKHelper".
# Note: Order matters below, due to inter-module dependencies.
include(joinpath("utils", "Serialization.jl"))
@include_export("NEPCore.jl", NEPCore)
@include_export("NEPTypes.jl", NEPTypes)
@include_export("LinSolvers.jl", LinSolvers)
include(joinpath("rk_helper", "RKHelper.jl"))
@include_export("NEPSolver.jl", NEPSolver)
@include_export("Gallery.jl", Gallery)

end
