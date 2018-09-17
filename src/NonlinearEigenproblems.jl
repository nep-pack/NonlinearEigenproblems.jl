module NonlinearEigenproblems

    include(joinpath("utils", "Serialization.jl"))
    include("NEPCore.jl")
    include("NEPTypes.jl")
    include("LinSolvers.jl")
    include("NEPSolver.jl")
    include("Gallery.jl")

    using Reexport
    @reexport using .NEPCore
    @reexport using .NEPTypes
    @reexport using .NEPSolver
    @reexport using .LinSolvers
    @reexport using .Gallery

end
