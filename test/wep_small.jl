# Run tests for the waveguide eigenvalue problem

using NonlinearEigenproblems
using Test
using LinearAlgebra

using GalleryWaveguide

import GalleryWaveguide.SchurMatVec

@bench @testset "WEP" begin

    @bench @testset "SPMF - NEP, and Preconditioner" begin
        @info "Compare SPFM and native format, and test Preconditioner"
        nx = 11
        nz = 7
        nep_spmf=nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = "TAUSCH", neptype = "SPMF")
        nep=nep_gallery(WEP, nx = nx, nz = nz, benchmark_problem = "TAUSCH", neptype = "WEP")
        λ = -1.3-0.31im
        v1 = compute_Mlincomb(nep_spmf, λ, ones(size(nep_spmf,1)))
        v2 = compute_Mlincomb(nep     , λ, ones(size(nep     ,1)))
        @test norm(v1-v2)/norm(v1) < 1e-14

        precond = wep_generate_preconditioner(nep, nz, λ)
        b1 = rand(ComplexF64, nx*nz)
        Schur_fun = SchurMatVec(nep, λ)
        b2 = ldiv!(precond, (Schur_fun*b1))
        @test norm(b1-b2)/norm(b1) < 1e-14
    end

    nep=nep_gallery(WEP, nx = 3*5*7+4, nz = 3*5*7, benchmark_problem = "JARLEBRING", neptype = "WEP")
    n=size(nep,1);
    λ0=-3-3.5im
    v0=ones(n); v0=v0/norm(v0);
    λref =-2.743228671961724 - 3.1439375599649972im  # Reference eigenvalue
    myerrmeasure = EigvalReferenceErrmeasure(nep,λref);


    @bench @testset "Linear solvers" begin

        @info "WEP default linsolver"
        λ,v = resinv(ComplexF64,nep,logger=displaylevel,λ=λ0,v=v0,
                     errmeasure=myerrmeasure,tol=1e-12,
                     linsolvercreator=GalleryWaveguide.WEPLinSolverCreator()
                     );
        @test  norm(compute_Mlincomb(nep,λ,v))/norm(v)  < 1e-10

        @info "WEP backslash linsolver"
        λ,v = resinv(ComplexF64,nep,logger=displaylevel,λ=λ0,v=v0,
                     errmeasure=myerrmeasure,tol=1e-12,
                     linsolvercreator=GalleryWaveguide.WEPLinSolverCreator(solver_type=:backslash)
                     );
        @test  norm(compute_Mlincomb(nep,λ,v))/norm(v)  < 1e-10

        @info "WEP gmres linsolver"
        precond = wep_generate_preconditioner(nep, 3*7, λ0)
        λ,v = resinv(ComplexF64,nep,logger=displaylevel,λ=λ0,v=v0,
                     errmeasure=myerrmeasure,tol=1e-12,
                     linsolvercreator=GalleryWaveguide.WEPLinSolverCreator(solver_type=:gmres,kwargs=((:Pl,precond),(:reltol,1e-7)))
                     );
        @test  norm(compute_Mlincomb(nep,λ,v))/norm(v)  < 1e-10
    end

    @bench @testset "NEP solvers" begin
        @info "WEP NEP solvers"
        λ,v = quasinewton(ComplexF64,nep,logger=displaylevel,λ=λ0,v=v0,
                          errmeasure=myerrmeasure,tol=1e-12,
                          linsolvercreator=GalleryWaveguide.WEPLinSolverCreator(solver_type=:factorized)
                          );
        @test  norm(compute_Mlincomb(nep,λ,v))/norm(v)  < 1e-10

        nev=3
        λ,v = iar(ComplexF64,nep,σ=λ0, logger=displaylevel,neigs=nev,maxit=100,v=v0,
                  tol=1e-8, linsolvercreator=GalleryWaveguide.WEPLinSolverCreator(solver_type=:factorized)
                  );
        @test minimum(abs.(λref .- λ)) < 1e-10
    end



end
