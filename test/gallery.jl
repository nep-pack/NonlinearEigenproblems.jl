# Run tests on nep_gallery (not tested elsewhere)

using NonlinearEigenproblems
using LinearAlgebra
using SparseArrays
using Test


@testset "Gallery" begin
    @testset "Gallery extra loading" begin
        # Test if it is possible to load some extra gallery
        @test (using GalleryWaveguide) == nothing
    end

    @testset "Basic random examples" begin
        @info "Testing dep0"
        nep = nep_gallery("dep0")
        nep = nep_gallery("dep0", 15)
        A = compute_Mder(nep, 3.0)
        @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery("dep0", 15, 8)
        @test_throws MethodError nep_gallery("dep0", 15, t=8)

        @info "Testing dep0_sparse"
        nep = nep_gallery("dep0_sparse")
        nep = nep_gallery("dep0_sparse", 20)
        nep = nep_gallery("dep0_sparse", 20, 0.25)
        A = compute_Mder(nep, 3.0)
        @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery("dep0_sparse", 20, 0.25, 8)
        @test_throws MethodError nep_gallery("dep0_sparse", 20, 0.25, n=8)

        @info "Testing dep0_tridiag"
        nep = nep_gallery("dep0_tridiag")
        nep = nep_gallery("dep0_tridiag", 15)
        A = compute_Mder(nep, 3.0)
        @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery("dep0", 15, 8)
        @test_throws MethodError nep_gallery("dep0", 15, t=8)

        @info "pep0"
        nep = nep_gallery("pep0")
        nep = nep_gallery("pep0", 15)
        A = compute_Mder(nep, 3.0)
        @test isa(nep, PEP) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery("pep0", 15, 8)
        @test_throws MethodError nep_gallery("pep0", 15, t=8)

        @info "pep0_sym"
        nep = nep_gallery("pep0_sym")
        nep = nep_gallery("pep0_sym", 15)
        A = compute_Mder(nep, 3.0)
        @test isa(nep, PEP) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery("pep0_sym", 15, 8)
        @test_throws MethodError nep_gallery("pep0_sym", 15, t=8)

        @info "Testing pep0_sparse"
        nep=nep_gallery("pep0_sparse")
        nep=nep_gallery("pep0_sparse", 15)
        nep=nep_gallery("pep0_sparse", 15, 0.12)
        A=compute_Mder(nep,3.0)
        @test issparse(A) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery("pep0_sparse", 15, 0.12, 8)
        @test_throws MethodError nep_gallery("pep0_sparse", 15, 0.12, t=8)

        @info "Testing pdde_stability"
        pep=nep_gallery("nlevp_native_pdde_stability",3);
        λ=-0.5064627366803098 - 0.00019049029618525343im; # Reference solution
        @test minimum(svdvals!(Matrix(compute_Mder(pep,λ)))) < eps()*100;

    end

    @info "Testing dep_symm_double"
    nep = nep_gallery("dep_symm_double")
    nep = nep_gallery("dep_symm_double", 15)
    A = compute_Mder(nep, 3.0)
    @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("dep_symm_double", 15, 8)
    @test_throws MethodError nep_gallery("dep_symm_double", 15, t=8)

    @info "Testing dep_double"
    nep = nep_gallery("dep_double")
    A = compute_Mder(nep, 3.0)
    @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("dep_double", 15)
    @test_throws MethodError nep_gallery("dep_double", t=15)

    @info "real_quadratic"
    nep = nep_gallery("real_quadratic")
    A = compute_Mder(nep, 3.0)
    @test isa(nep, PEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("real_quadratic", 15)
    @test_throws MethodError nep_gallery("real_quadratic", t=15)


    for i = [0, 1]
        nep_name = "qdep" * string(i)
        @info nep_name
        nep = nep_gallery(nep_name)
        A = compute_Mder(nep, 3.0)
        @test isa(nep, SPMF_NEP) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery(nep_name, 15)
        @test_throws MethodError nep_gallery(nep_name, t=15)
    end

    @info "Testing sine"
    nep=nep_gallery("sine")
    tol=1e-10
    λ0=λ=-1.4;
    v0=normalize!(compute_Mder(nep,λ0)\ones(size(nep,1)))
    λ,v=quasinewton(Float64,nep,λ=-1.4566,v=v0,tol=tol,
                    errmeasure=ResidualErrmeasure(nep),
                    logger=displaylevel,armijo_factor=0.5,armijo_max=3)
    normalize!(v)
    @test norm(compute_Mlincomb(nep,λ,v))<sqrt(tol)
    @test isa(nep,SPMFSumNEP)
    @test_throws MethodError nep_gallery("sine", 15)
    @test_throws MethodError nep_gallery("sine", t=15)


    @info "Testing dep_symm_double"
    nep=nep_gallery("dep_symm_double")
    A=compute_Mder(nep,3.0)
    @test norm(A-A')==0

    @info "Testing qep_fixed_eig"
    nep=nep_gallery("qep_fixed_eig")
    nep=nep_gallery("qep_fixed_eig",3)
    nep=nep_gallery("qep_fixed_eig",3,1:6)
    A=compute_Mder(nep,3.0)
    λstar=5.0 # this problem is constructed such that this is an eigenvalue
    @test minimum(svdvals(compute_Mder(nep,λstar)))<100*eps()
    @test_throws MethodError nep_gallery("qep_fixed_eig",3,1:6, 8)
    @test_throws MethodError nep_gallery("qep_fixed_eig",3,1:6, t=8)

    @info "Testing neuron0"
    nep=nep_gallery("neuron0")
    A=compute_Mder(nep,3.0)
    @test_throws MethodError nep_gallery("neuron0", 8)
    @test_throws MethodError nep_gallery("neuron0", t=8)

    @info "Testing schrodinger_movebc"
    n=5; λ=-3.0;
    nep=nep_gallery("schrodinger_movebc",n)
    A=compute_Mder(nep,λ)
    h=1/(n-1);
    @test A[1,1]== -2/(h^2)-λ-1 # Check (1,1)-point in discretization

    @info "Testing beam"
    nep=nep_gallery("beam")
    nep=nep_gallery("beam", 15)
    A=compute_Mder(nep,3.0)
    @test_throws MethodError nep_gallery("beam", 15, 8)
    @test_throws MethodError nep_gallery("beam", 15, t=8)


    @info "Testing hadeler"
    nep=nep_gallery("nlevp_native_hadeler")
    A=compute_Mder(nep,3.0);
    @test isreal(A)


    @info "Testing loaded_string"
    # Compare with hard-coded value
    λ=3.0+4im;
    Zref=[9.600000000000000 - 0.533333333333333im -5.100000000000000 - 0.133333333333333im  0.000000000000000 + 0.000000000000000im  0.000000000000000 + 0.000000000000000im  0.000000000000000 + 0.000000000000000im
 -5.100000000000000 - 0.133333333333333im  9.600000000000000 - 0.533333333333333im -5.100000000000000 - 0.133333333333333im  0.000000000000000 + 0.000000000000000im  0.000000000000000 + 0.000000000000000im
  0.000000000000000 + 0.000000000000000im -5.100000000000000 - 0.133333333333333im  9.600000000000000 - 0.533333333333333im -5.100000000000000 - 0.133333333333333im  0.000000000000000 + 0.000000000000000im
  0.000000000000000 + 0.000000000000000im  0.000000000000000 + 0.000000000000000im -5.100000000000000 - 0.133333333333333im  9.600000000000000 - 0.533333333333333im -5.100000000000000 - 0.133333333333333im
       0.000000000000000 + 0.000000000000000im  0.000000000000000 + 0.000000000000000im  0.000000000000000 + 0.000000000000000im -5.100000000000000 - 0.133333333333333im  9.137811900191938 - 0.880870121561100im ];
    nep=nep_gallery("nlevp_native_loaded_string",5,4,5);
    Z=compute_Mder(nep,λ);
    @test norm(Z-Zref)<1e-10




    @info "testing bem_fichera"
    nep=nep_gallery("bem_fichera",1);

    λ_ref=8.790558462139456 - 0.010815457827738698im
    M=compute_Mder(nep,λ_ref);
    @test rank(M)<size(M,1)
    @onlybench @testset "bem_fichera + IAR" begin
       v=ones(size(nep,1));
       (λ,vv)=iar(nep,σ=9,v=v,logger=displaylevel,neigs=1,maxit=50,tol=1e-6) # normally takes 30 iterations
       @test norm(compute_Mlincomb(nep,λ[1],vv[:,1]))<sqrt(eps())*100;
    end
    @info "non-existing example"
    @test_throws ErrorException nep_gallery("non-existing example")

end
