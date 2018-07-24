# Test of companion linearization

if !isdefined(:global_running_all_tests) || global_running_all_tests != true
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using Base.Test
end


comaniontest = @testset "Companion Linearization" begin

#####################
# Dense matrix test #
#####################
pep = nep_gallery("pep0");
deg = size(get_fv(pep),1) -1;
n = size(pep,1);
tolerance = 1e-14;

λa,xa = newton(pep, maxit=30, displaylevel=0, tol = tolerance);
Dc,Vc = polyeig(pep, DefaultEigSolver);

ind = 1;
for i = 1:(deg*n)
    if(abs(λa-Dc[i])/abs(λa) < tolerance*10)
        ind = i;
        break
    end
end
@test (abs(λa-Dc[ind])/abs(λa) < tolerance*10)


######################
# Sparse matrix test #
######################
pep = nep_gallery("pep0_sparse_003");
deg = size(get_fv(pep),1) -1;
n = size(pep,1);
tolerance = 1e-14;

λa,xa =newton(pep, maxit=30, λ=-0.75, v=ones(n), displaylevel=0, tol = tolerance);
Dc,Vc = polyeig(pep,DefaultEigSolver);

ind = 1;
for i = 1:(deg*n)
    if(abs(λa-Dc[i])/abs(λa) < tolerance*10)
        ind = i;
        break
    end
end
@test (abs(λa-Dc[ind])/abs(λa) < tolerance*10)


########################
# BigFloat matrix test #
########################
A0 = (Array{BigFloat,2}([1 2; 3 4]));
A1 = (Array{BigFloat,2}([1 44; 3 4]));
A2 = (Array{BigFloat,2}([1 44; -3 100]));
pep = PEP([A0,A1,A2])

# Power method for testing (since eig does not work for BigFloats)
E,A = companion(pep);
z=ones(BigFloat,size(A,1));
local evp::BigFloat
local evp_old::BigFloat
local d::BigFloat
tolerance = big"1e-50"
d=Inf; evp=Inf
k=0;
while abs(d) > tolerance
    k=k+1
    z=z/norm(z);
    z2=E\A*z
    evp_old=evp
    evp=dot(z,z2)
    d=evp-evp_old
    # println("Iteration ",k,", λ=",Float64(evp), " Δ = ",Float64(d))
    z=z2
end

println("Solving same problem with resinv")
λ,v=resinv(BigFloat, pep, λ = (BigFloat(evp)+0.1), v = z[1:size(pep,1)], displaylevel=0, tol = tolerance)

@test (abs(λ-evp)/abs(λ) < tolerance*10)


end
