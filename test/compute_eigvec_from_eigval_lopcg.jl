if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))
    push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
    push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using IterativeSolvers
    using Base.Test
end

# explicit import needed for overloading functions from packages
import NEPCore.compute_Mlincomb

@testset "dep0_sparse" begin
    nep=nep_gallery("dep0_sparse",100);nept=DEP([nep.A[1]',nep.A[2]'],nep.tauv);
    λ,Q,err = iar(nep,maxit=100,Neig=2,σ=1.0,γ=1,displaylevel=0,check_error_every=1);
    v=compute_eigvec_from_eigval_lopcg(nep,nept,λ[1]);
    errormeasure=default_errmeasure(nep);
    @test errormeasure(λ[1],v)<1e-5;
end

@testset "pep0" begin
    nep=nep_gallery("pep0",100);
    nept=PEP([nep.A[1]',nep.A[2]',nep.A[3]'])
    λ,Q,err = iar(nep,maxit=100,Neig=2,σ=1.0,γ=1,displaylevel=0,check_error_every=1);
    v=compute_eigvec_from_eigval_lopcg(nep,nept,λ[1]);
    errormeasure=default_errmeasure(nep);
    @test errormeasure(λ[1],v)<1e-5;
end
