# Run tests for the waveguide eigenvalue problem

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
#push!(LOAD_PATH, pwd()*"/src/gallery_extra")
#push!(LOAD_PATH, pwd()*"/src/gallery_extra/waveguide")	

using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers
using Gallery.Waveguide

using Base.Test
import NEPCore.compute_Mder;


nep_spmf=nep_gallery("waveguide", nx = 3*5*7, nz = 3*5*7, benchmark_problem = "JARLEBRING", discretization = "FD", neptype = "SPMF")

AA=nep_spmf.A;
AAt=Array{SparseMatrixCSC,1}(size(AA,1));
for i=1:size(AA,1)
    AAt[i]=AA[i]';
end

nept_spmf=SPMF_NEP(AAt,nep_spmf.fi)




nep=nep_gallery("waveguide", nx = 3*5*7, nz = 3*5*7, benchmark_problem = "JARLEBRING", discretization = "FD", neptype = "WEP")
# Waveguide compute_Mder(WEP_FD) (or direct solve) not implemented.
# Workaround: use the spmf-version, but only for the linear solves.
nep_spmf=nep_gallery("waveguide", nx = 3*5*7, nz = 3*5*7, benchmark_problem = "JARLEBRING", discretization = "FD", neptype = "SPMF")
function compute_Mder(nep::WEP_FD,λ::Number,i::Integer)
    return compute_Mder(nep_spmf,λ,i);
end
λstar=-2.690050173308845 - 3.1436003386330347im  # An exact eigenvalue
n=size(nep,1);

@testset "WEP" begin
    λ0=-3-3.5im
    v0=ones(n); v0=v0/norm(v0);
    #
    n0=norm(compute_Mlincomb(nep,λ0,v0))
    #myerrmeasure=(λ,v) -> (norm(compute_Mlincomb(nep,λ,v))/(n0*norm(v)))
    myerrmeasure=(λ,v) -> abs(λ-λstar) # Use eigenvalue error as errmeasure

    λ,v=resinv(Complex128,nep,displaylevel=1,λ=λ0,v=v0,
               errmeasure=myerrmeasure,tol=1e-12)


    @test  norm(compute_Mlincomb(nep,λ,v))/norm(v)  < 1e-10
    λ,v=quasinewton(Complex128,nep,displaylevel=1,λ=λ0,v=v0,
                    errmeasure=myerrmeasure,tol=1e-12)

    @test  norm(compute_Mlincomb(nep,λ,v))/norm(v)  < 1e-10


    nev=3
   
    λ,v=@time iar(Complex128,nep,σ=λ0, displaylevel=1,Neig=nev,maxit=100,v=v0,
                  tol=1e-8);

    @test minimum(abs.(λstar-λ))<1e-10

    #λ,v=@time tiar(nep,σ=λ0, displaylevel=1,Neig=nev,maxit=100,v=v0,
    #               tol=1e-8);
    
end


