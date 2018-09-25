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


@bench @testset "Gallery" begin
    @info "Testing sine"
    nep=nep_gallery("sine")
    tol=1e-10;
    λ,v=quasinewton(Float64,nep,λ=-4.2,v=ones(size(nep,1)),tol=tol,
                    displaylevel=displaylevel,armijo_factor=0.5,armijo_max=20)
    normalize!(v);
    @test norm(compute_Mlincomb(nep,λ,v))<tol*100

    @info "Testing dep_symm_double"
    nep=nep_gallery("dep_symm_double");
    A=compute_Mder(nep,3.0);
    @test norm(A-A')==0;


    @info "Testing pep0_sparse_003"
    nep=nep_gallery("pep0_sparse_003");;
    A=compute_Mder(nep,3.0);
    @test issparse(A); # Not so sophisticated, but improves code coverage


    @info "Testing qep_fixed_eig"
    nep=nep_gallery("qep_fixed_eig",3,1:6);;
    A=compute_Mder(nep,3.0);
    λstar=5.0; # this problem is constructed such that this is an eigenvalue
    @test minimum(svdvals(compute_Mder(nep,λstar)))<100*eps()

end
