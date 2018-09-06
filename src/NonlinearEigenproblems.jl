module NonlinearEigenproblems

include(joinpath("utils", "Serialization.jl"))
include("NEPCore.jl")
include("NEPTypes.jl")
include("LinSolvers.jl")
include("NEPSolver.jl")
include("Gallery.jl")

end
