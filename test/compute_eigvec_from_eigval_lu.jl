if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

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

nep_test_problems=["pep0_sparse_003","dep0","pep0"]

eigeinvector_extraction_small=@testset "Eigenvector extraction (small scale)" begin
    @testset "Test problem: $nep_test_problem" for nep_test_problem in nep_test_problems
    nep=nep_gallery(nep_test_problem)
    compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
    errormeasure=default_errmeasure(nep);
    λ,Q,err = iar(nep,maxit=50,Neig=5,σ=2.0,γ=3);
        @testset "default_linsolvercreator" begin
            v=compute_eigvec_from_eigval_lu(nep,λ[1],default_linsolvercreator);
            @test errormeasure(λ[1],v)<eps()*10000;
        end
        @testset "backslash_linsolvercreator" begin
            M0inv = backslash_linsolvercreator(nep,λ[1])
            v=compute_eigvec_from_eigval_lu(nep,λ[1],default_linsolvercreator);
            @test errormeasure(λ[1],v)<eps()*10000;
        end
        @testset "Passing a linsolvercreator as argument" begin
            M0inv = default_linsolvercreator(nep,λ[1])
            v=compute_eigvec_from_eigval_lu(nep,λ[1],(nep, σ) -> M0inv);
            @test errormeasure(λ[1],v)<eps()*10000;
        end
    end
end


eigeinvector_extraction_large=@testset "Eigenvector extraction (medium/large scale)" begin
    @testset "Test problem: $nep_test_problem" for nep_test_problem in nep_test_problems
    nep=nep_gallery(nep_test_problem,500)
    compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
    errormeasure=default_errmeasure(nep);
    λ,Q,err = iar(nep,maxit=100,Neig=5);
        @testset "default_linsolvercreator" begin
            v=compute_eigvec_from_eigval_lu(nep,λ[1],default_linsolvercreator);
            @test errormeasure(λ[1],v)<eps()*10000;
        end
        @testset "backslash_linsolvercreator" begin
            M0inv = backslash_linsolvercreator(nep,λ[1])
            v=compute_eigvec_from_eigval_lu(nep,λ[1],default_linsolvercreator);
            @test errormeasure(λ[1],v)<eps()*10000;
        end
        @testset "Passing a linsolvercreator as argument" begin
            M0inv = default_linsolvercreator(nep,λ[1])
            v=compute_eigvec_from_eigval_lu(nep,λ[1],(nep, σ) -> M0inv);
            @test errormeasure(λ[1],v)<eps()*10000;
        end
    end
end
