# Run tests on block Newton

# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_running_all_tests) || global_running_all_tests != true
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using Base.Test
end

@testset "sgiter" begin

    nep = nep_gallery("real_quadratic")
#    nep = nep_gallery("dep_distributed")
#    nep = nep_gallery("pep0_sym")
    λ,v = sgiter(Float64,
                 nep,
                 1,
                 λ_min = -10,
                 λ_max = 0,
                 λ = -10,
                 displaylevel = 1,
                 maxit = 100)
   @test norm(compute_Mlincomb(nep,λ,v)) < eps(Float64)*100


   λ,v = sgiter(Float64,
                nep,
                2,
                displaylevel = 1,
                tol = 1e-9,
                maxit = 100)
  @test norm(compute_Mlincomb(nep,λ,v)) < 1e-9

end
