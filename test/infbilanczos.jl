# Test for infinite Bi-Lanczos

workspace()
push!(LOAD_PATH, pwd()*"/src")	
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using MAT
using Base.Test
using LinSolvers


# explicit import needed for overloading
# functions from packages
# import NEPCore.compute_Mlincomb

#println("Load dep0")
#nep=nep_gallery("dep0",100)
##nep=nep_gallery("pep0");
#println("Transposing dep0");
## Construct the transposed NEP
#nept=DEP([nep.A[1]',nep.A[2]'],nep.tauv)
#
#


file=matopen("matlab/infbilanczos/qdep_infbilanczos.mat");
A0=read(file,"A0");
A1=read(file,"A1");
close(file)

tau=1;
quadfun= S -> S^2;
constfun= S -> eye(S);
expfun= S -> expm(-tau*S);


AA=[-speye(A0),A0,A1]

fi=[quadfun,constfun,expfun]
nep=SPMF_NEP(AA,fi)

AAt=[-speye(A0),A0',A1']

n=size(nep,1);

nept=SPMF_NEP(AAt,[quadfun,constfun,expfun])


@testset "Infbilanczos σ=0" begin
    m=30;
    λ,Q,T = infbilanczos(Float64,nep,nept,maxit=m,Neig=m,σ=0,
                         v=ones(n),u=ones(n));

    # Produced with a different implementation
    Tstar=[  -1.665117675679600   5.780562035399026                   0                   0
             5.780562035399026  11.562308485001218 -18.839546184493731                   0
             0  18.839546184493734 -15.213756300995186   9.788512505128466                                
             0                   0   9.788512505128464  -0.120825360586847                                ] 
    n0=size(Tstar,1);
    @test norm(Tstar-T[1:n0,1:n0])<1e-10

    println("Post-processing: computing eigenvectors")
    thiserr=ones(m)*NaN
    for i=1:length(λ)
        v=compute_eigvec_from_eigval_lu(nep,λ[i],default_linsolvercreator);
        thiserr[i]=norm(compute_Mlincomb(nep,λ[i],v));
    end
    @test length(find(thiserr .< 1e-7)) == 9
end


@testset "Infbilanczos σ=$x" for x in (1.0, 1.0+0.2im)
    m=20;
    λ,Q,T = infbilanczos(nep,nept,maxit=m,Neig=m,σ=x,
                         v=ones(n),u=ones(n));
    thiserr=ones(m)*NaN
    for i=1:length(λ)
        v=compute_eigvec_from_eigval_lu(nep,λ[i],default_linsolvercreator);
        thiserr[i]=norm(compute_Mlincomb(nep,λ[i],v));
    end
    @test length(find(thiserr .< 1e-7)) >= 3
end

    




#
##for i=1:m
## semilogy(1:m, err[1:m,i], color="red", linewidth=2.0, linestyle="--")
##end
##
#errormeasure=default_errmeasure(nep);
#for i=1:length(λ)
# println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
#end
#
