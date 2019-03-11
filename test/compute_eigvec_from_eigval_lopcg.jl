using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Random
using Test

# explicit import needed for overloading functions from packages
import NonlinearEigenproblems.NEPCore.compute_Mlincomb

@testset "compute eigvec lopcg" begin
    # Arbitrary vectors. We avoid random in unit tests
    v0=Vector{Float64}(cos.(3.2 .+ (1:100))); normalize!(v0);
    v1=Vector{Float64}(sin.(1:100)); normalize!(v0);
    @bench @testset "dep0_sparse" begin
        nep = nep_gallery("dep0_sparse", 100);
        nept = DEP([copy(nep.A[1]'), copy(nep.A[2]')], nep.tauv);
        λ,Q,err = iar(nep,maxit=100,Neig=2,σ=1.0,γ=1,displaylevel=0,check_error_every=1,v=v0);
        v=compute_eigvec_from_eigval_lopcg(nep,nept,λ[1],x=v1);
        ermdata=init_errmeasure(ResidualErrmeasure,nep);
        @test estimate_error(ermdata,λ[1],v)<1e-5;
    end

    @bench @testset "pep0" begin
        nep = nep_gallery("pep0", 100);
        nept = PEP([copy(nep.A[1]'), copy(nep.A[2]'), copy(nep.A[3]')])
        λ,Q,err = iar(nep,maxit=100,Neig=2,σ=1.0,γ=1,displaylevel=0,check_error_every=1,v=v0);
        v=compute_eigvec_from_eigval_lopcg(nep,nept,λ[1],x=v1);
        ermdata=init_errmeasure(ResidualErrmeasure,nep);
        @test estimate_error(ermdata,λ[1],v)<1e-5;
    end

end
