using NonlinearEigenproblems
using Test
using LinearAlgebra

@bench @testset "sgiter" begin

    nep = nep_gallery("real_quadratic")
#    nep = nep_gallery("dep_distributed")
#    nep = nep_gallery("pep0_sym")
    λ,v = sgiter(Float64,
                 nep,
                 1,
                 λ_min = -10,
                 λ_max = 0,
                 λ = -10,
                 logger = displaylevel,
                 maxit = 100, tol=1e-12)
   @test norm(compute_Mlincomb(nep,λ,v)) < 1e-9


   λ,v = sgiter(Float64,
                nep,
                2,
                logger = displaylevel,
                tol = 1e-9,
                maxit = 100,
                errmeasure = ResidualErrmeasure(nep))
  @test norm(compute_Mlincomb(nep,λ,v)) < 1e-9

end
