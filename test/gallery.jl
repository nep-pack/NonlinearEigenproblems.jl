# Run tests on nep_gallery (not tested elsewhere)

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems
using LinearAlgebra
using SparseArrays
using Test

@testset "Gallery extra loading" begin
    # Test if it is possible to load some extra gallery
    @test (using GalleryWaveguide) == nothing
end


@testset "Gallery" begin
    @info "Testing dep0"
    nep = nep_gallery("dep0")
    nep = nep_gallery("dep0", 15)
    A = compute_Mder(nep, 3.0)
    @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("dep0", 15, 8)
    @test_throws ErrorException nep_gallery("dep0", 15, t=8)

    @info "Testing dep0_sparse"
    nep = nep_gallery("dep0_sparse")
    nep = nep_gallery("dep0_sparse", 20)
    nep = nep_gallery("dep0_sparse", 20, 0.25)
    A = compute_Mder(nep, 3.0)
    @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("dep0_sparse", 20, 0.25, 8)
    @test_throws ErrorException nep_gallery("dep0_sparse", 20, 0.25, n=8)

    @info "Testing dep0_tridiag"
    nep = nep_gallery("dep0_tridiag")
    nep = nep_gallery("dep0_tridiag", 15)
    A = compute_Mder(nep, 3.0)
    @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("dep0", 15, 8)
    @test_throws ErrorException nep_gallery("dep0", 15, t=8)

    @info "Testing dep_symm_double"
    nep = nep_gallery("dep_symm_double")
    nep = nep_gallery("dep_symm_double", 15)
    A = compute_Mder(nep, 3.0)
    @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("dep_symm_double", 15, 8)
    @test_throws ErrorException nep_gallery("dep_symm_double", 15, t=8)

    @info "Testing dep_double"
    nep = nep_gallery("dep_double")
    A = compute_Mder(nep, 3.0)
    @test isa(nep, DEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("dep_double", 15)
    @test_throws ErrorException nep_gallery("dep_double", t=15)

    @info "pep0"
    nep = nep_gallery("pep0")
    nep = nep_gallery("pep0", 15)
    A = compute_Mder(nep, 3.0)
    @test isa(nep, PEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("pep0", 15, 8)
    @test_throws ErrorException nep_gallery("pep0", 15, t=8)

    @info "pep0_sym"
    nep = nep_gallery("pep0_sym")
    nep = nep_gallery("pep0_sym", 15)
    A = compute_Mder(nep, 3.0)
    @test isa(nep, PEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("pep0_sym", 15, 8)
    @test_throws ErrorException nep_gallery("pep0_sym", 15, t=8)

    @info "real_quadratic"
    nep = nep_gallery("real_quadratic")
    A = compute_Mder(nep, 3.0)
    @test isa(nep, PEP) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("real_quadratic", 15)
    @test_throws ErrorException nep_gallery("real_quadratic", t=15)


    for i = [0, 1]
        nep_name = "qdep" * string(i)
        @info nep_name
        nep = nep_gallery(nep_name)
        A = compute_Mder(nep, 3.0)
        @test isa(nep, SPMF_NEP) # Not so sophisticated, but improves code coverage
        @test_throws MethodError nep_gallery(nep_name, 15)
        @test_throws ErrorException nep_gallery(nep_name, t=15)
    end

    @info "Testing sine"
    nep=nep_gallery("sine")
    tol=1e-10
    λ,v=quasinewton(Float64,nep,λ=-4.2,v=ones(size(nep,1)),tol=tol,
                    displaylevel=displaylevel,armijo_factor=0.5,armijo_max=20)
    normalize!(v)
    @test norm(compute_Mlincomb(nep,λ,v))<tol*100

    @info "Testing dep_symm_double"
    nep=nep_gallery("dep_symm_double")
    A=compute_Mder(nep,3.0)
    @test norm(A-A')==0


    @info "Testing pep0_sparse"
    nep=nep_gallery("pep0_sparse")
    nep=nep_gallery("pep0_sparse", 15)
    nep=nep_gallery("pep0_sparse", 15, 0.12)
    A=compute_Mder(nep,3.0)
    @test issparse(A) # Not so sophisticated, but improves code coverage
    @test_throws MethodError nep_gallery("pep0_sparse", 15, 0.12, 8)
    @test_throws ErrorException nep_gallery("pep0_sparse", 15, 0.12, t=8)


    @info "Testing qep_fixed_eig"
    nep=nep_gallery("qep_fixed_eig")
    nep=nep_gallery("qep_fixed_eig",3)
    nep=nep_gallery("qep_fixed_eig",3,1:6)
    A=compute_Mder(nep,3.0)
    λstar=5.0 # this problem is constructed such that this is an eigenvalue
    @test minimum(svdvals(compute_Mder(nep,λstar)))<100*eps()
    @test_throws MethodError nep_gallery("qep_fixed_eig",3,1:6, 8)
    @test_throws ErrorException nep_gallery("qep_fixed_eig",3,1:6, t=8)

    @info "Testing neuron0"
    nep=nep_gallery("neuron0")
    A=compute_Mder(nep,3.0)
    @test_throws MethodError nep_gallery("neuron0", 8)
    @test_throws ErrorException nep_gallery("neuron0", t=8)

end
