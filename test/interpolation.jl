using NonlinearEigenproblems
using Test
using LinearAlgebra


@testset "Interpolation" begin


#########################################################################################
@bench @testset "Random dep" begin
nep = nep_gallery("dep0")
λ1,x1 = newton(nep,logger=displaylevel,maxit=40, λ=-0.75, v=ones(size(nep,1)));
@test compute_resnorm(nep,λ1,x1) < eps()*100

# Newton on interpolated dep (interpolating pep of degree 2)
intpoints = [λ1-1, λ1, λ1+1]
pep = interpolate(nep, intpoints)
λ2,x2 = newton(pep,logger=displaylevel,maxit=40, λ=-0.75, v=ones(size(nep,1)));
@test compute_resnorm(pep,λ2,x2) < eps()*100
@test abs(λ1-λ2)/abs(λ1) < eps()*1000

# Newton on interpolated dep (interpolating pep of degree 8)
intpoints = [λ1-5, λ1-1, λ1, λ1+5, λ1+1, λ1+5im, λ1+1im, λ1-5im, λ1-1im]
pep = interpolate(nep, intpoints)
λ2,x2 = newton(pep,logger=displaylevel,maxit=40, λ=-0.75, v=ones(size(nep,1)));
@test compute_resnorm(pep,λ2,x2) < eps()*100
@test abs(λ1-λ2)/abs(λ1) < eps()*1000
end

#########################################################################################
@bench @testset "Random sparse dep" begin
nep = nep_gallery("dep0_sparse", 30)

λ1,x1 = newton(nep,logger=displaylevel,maxit=40, λ=0.2, v=ones(size(nep,1)),armijo_factor=0.5);
@test compute_resnorm(nep,λ1,x1) < eps()*100

# Newton on interpolated sparse dep
intpoints = [λ1-.5, λ1, λ1+.5]
pep = interpolate(nep, intpoints)
λ2,x2 = newton(pep,logger=displaylevel,maxit=40, λ=λ1+0.1, v=x1+0.1*ones(size(nep,1)));
@test compute_resnorm(pep,λ2,x2) < eps()*1000000
@test abs(λ1-λ2)/abs(λ1) < eps()*1000

end

#########################################################################################
@bench @testset "Random pep (original pep of degree 2)" begin
nep = nep_gallery("pep0")
λ1,x1 = newton(nep,logger=displaylevel,maxit=40, λ=1, v=ones(size(nep,1)));
@test compute_resnorm(nep,λ1,x1) < eps()*100


#Newton on interpolated pep (interpolation of degree 2)
intpoints = [λ1-1, λ1, λ1+1]
pep = interpolate(nep, intpoints)
λ2,x2 = newton(pep,logger=displaylevel,maxit=40, λ=1, v=ones(size(nep,1)));
@test compute_resnorm(pep,λ2,x2) < eps()*100
@test abs(λ1-λ2)/abs(λ1) < eps()*100

@testset "Coefficient matrix differences (monomes) deg 2" begin
for i = 1:3
    @test opnorm(nep.A[i]-pep.A[i])/opnorm(nep.A[i]) < eps()*100
end
end


#Newton on interpolated pep (interpolation of degree 4) in real arithmetics
intpoints = [λ1-3, λ1-1, λ1, λ1+1, λ1+3]

pep = interpolate(Float64, nep, intpoints)
λ2,x2 = newton(pep,logger=displaylevel,maxit=40, λ=1, v=ones(size(nep,1)));
@test compute_resnorm(pep,λ2,x2) < eps()*100
@test abs(λ1-λ2)/abs(λ1) < eps()*100

@testset "Coefficient matrix differences (monomes) deg 4" begin
for i = 1:3
    @test opnorm(nep.A[i]-pep.A[i])/opnorm(nep.A[i]) < eps()*100
end
for i = 4:5
    @test opnorm(pep.A[i]) < eps()*100
end
end

end

#########################################################################################
@bench @testset "Random sparse pep (original pep of degree 2)" begin
nep=nep_gallery("pep0_sparse")
λ1,x1 =newton(nep,logger=displaylevel,maxit=40, λ=-0.75, v=ones(size(nep,1)));
@test compute_resnorm(nep,λ1,x1) < eps()*1000


#Newton on interpolated sparse pep
intpoints = [λ1-1, λ1, λ1+1.5]
pep = interpolate(nep, intpoints)
λ2,x2 =newton(pep,logger=displaylevel,maxit=40, λ=-0.75, v=ones(size(nep,1)));
@test compute_resnorm(pep,λ2,x2) < eps()*1000
@test abs(λ1-λ2)/abs(λ1) < eps()*100

@testset "Coefficient matrix differences (monomes) deg 2" begin
for i = 1:3
 @test opnorm(nep.A[i]-pep.A[i],Inf)/opnorm(nep.A[i],Inf) < eps()*100
end
end

end

end
