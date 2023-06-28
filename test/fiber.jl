# Run tests on the fiber problem in NLEVP (bessel function nonlinearity)

using NonlinearEigenproblems
using Test
using LinearAlgebra

push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryNLEVP

@testset "NLEVP fiber" begin
    nep_org=nep_gallery(NLEVP_NEP,"fiber");
    n=size(nep_org,1);

    # An exact eigenvalue according (reported in NLEVP collection)
    sol_val= 7.139494306065948e-07;

    @bench @testset "basics" begin
        # Checking Mder
        λ=3;
        δ=0.0001;
        Ap=compute_Mder(nep_org,λ+δ)
        Am=compute_Mder(nep_org,λ-δ)
        Mp_approx=(Ap-Am)/(2δ);
        Mp=compute_Mder(nep_org,λ,1)
        @test opnorm(Mp-Mp_approx,1)<δ^2*100


        # Checking Mlincomb
        λ=2;
        x=randn(n)
        z_plus=compute_Mlincomb(nep_org,λ+δ,x)
        z_minus=compute_Mlincomb(nep_org,λ-δ,x)
        zp=compute_Mlincomb(nep_org,λ,x,[1],1)
        zp-(z_plus-z_minus)/(2δ)
        @test norm(zp-(z_plus-z_minus)/(2δ))<(δ^2*100)


        @info "Running Newton"
        (λstar,v)=newton(Float64,nep_org,λ=7e-7,v=ones(n),logger=displaylevel,
                         tol=1e-14)

        @test abs(sol_val-λstar)/abs(λ) < 1e-10
        if imag(λstar) != 0
            @warn "Newton switches to complex although it should be real: $λstar"
        end

        @info "Running quasi-newton w armijo (with eigval error as termination)"
        function myerrmeasure(λt,v)
            return abs((λt-λstar)/λstar)
        end

        tol=1e-11
        (λ,v)=quasinewton(nep_org,λ=7.1e-7,v=ones(n),
                          logger=displaylevel,errmeasure=myerrmeasure,
                          tol=tol,armijo_factor=0.5,armijo_max=10)

        @test abs(sol_val-λ)/abs(λ) < tol
        if imag(λ) != 0
            @warn "Quasinewton switches to complex although it should be real: $λ"
        end

        @info "Running resinv"
        tol=1e-8
        (λ,v)=resinv(nep_org,λ=7e-7,v=ones(n),
                     logger=displaylevel,errmeasure=myerrmeasure,
                     tol=tol)
        @test abs(sol_val-λ)/abs(λ) < tol
        if imag(λ) != 0
            @warn "resinv switches to complex although it should be real: $λ"
        end

        @info "Running MSLP"

        (λ,v)=mslp(Float64,nep_org,λ=7e-7,
                   logger=displaylevel,errmeasure=myerrmeasure,
                   tol=tol)
        @test abs(sol_val-λ)/abs(λ) < tol
        if imag(λ) != 0
            @warn "mslp switches to complex although it should be real: $λ"
        end
    end



    @bench @testset "PEP interpolation" begin
        # Create a PEP Which interpolates in a couple of points

        intpoints=[1e-8,1e-7,7e-7,9e-7]
        @info "Interpolating: $intpoints"
        pep=interpolate(nep_org,intpoints);
        @info "Running IAR"
        λ,v=iar(pep,σ=7e-7,maxit=100,logger=displaylevel,neigs=2)
        minerr1=minimum(abs.(sol_val.-λ)) ./ abs(sol_val)
        @info "Error: $minerr1"
        @test minerr1<1e-4


        # Create a new PEP which interpolates in the solution computed
        # in the previous iteration
        intpoints=vcat(intpoints,real(λ[1]))
        @info "Interpolating: $intpoints"
        pep=interpolate(nep_org,intpoints);
        @info "Running IAR"
        λ,v=iar(pep,σ=7e-7,maxit=100,logger=displaylevel,neigs=2)
        minerr2=minimum(abs.(sol_val.-λ)) ./ abs(sol_val)
        @info "Error: $minerr2"
        @test minerr2<minerr1
    end

end
