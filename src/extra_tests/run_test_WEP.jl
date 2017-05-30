#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers


### DEBUGG ###
#import NEPCore.compute_Mder_from_MM
#import NEPCore.jordan_matrix
#export compute_Mder_from_MM
#function compute_Mder_from_MM(nep::SPMF_NEP,λ::Number,i::Integer=0)
#println("HERE: compute_Mder_from_MM(SPMF_NEP)")
#    J=sparse(jordan_matrix(typeof(λ),i+1,λ).')
#    n=size(nep,1);
#    S=kron(J,speye(n))
#    V=factorial(i)*kron(speye(1,i+1)[:,end:-1:1],speye(n))
#    @time W=compute_MM(nep,S,V)
#println("BYE!")
#    return W[1:n,1:n]
#end
#
#
#import NEPCore.compute_MM
#function compute_MM(nep::SPMF_NEP,S,V)
#        if (issparse(V))
#            if (size(V)==size(nep))
#                # Initialize with zero sparse matrix which 
#                # has sparsity pattern already consistent
#                # with sparsity pattern of M() for optimization
#                Z=copy(nep.Zero) 
#            else
#                Z=spzeros(eltype(V),size(V,1),size(V,2))
#            end
#        else
#            Z=zeros(eltype(V),size(V))
#        end
#        # Sum together all the terms in the SPMF:
#        @time for i=1:length(nep.A)
#            ## Compute Fi=f_i(S) in an optimized way
#            if (isdiag(S)) # optimize if S is diagonal
#                Sd=diag(S);
#                if (norm(Sd-Sd[1])==0) # Optimize further if S is a
#                                       # multiple of identity
#                    Fid=nep.fi[i](reshape([Sd[1]],1,1))[1]*ones(size(Sd,1))
#                else  # Diagonal but not constant
#                    Fid=zeros(Complex128,size(S,1))
#                    for j=1:size(S,1)
#                        Fid[j]=nep.fi[i](reshape([Sd[j]],1,1))[1]
#                    end                    
#                end
#                Fi=spdiagm(Fid);
#            else  # Otherwise just compute the matrix function operation
#                @time Fi=nep.fi[i](S)
#println("Just out")
#            end
#            ## Sum it up
#            Z=Z+nep.A[i]*(V*Fi);
#        end
#println("MM loop done\n\n\n\n\n")
#        return Z
#    end
### END ###






println("===========================")
println("||   This is WEP-test    ||")
println("===========================")


nz = 105;
nx = nz + 4;
delta = 0.1;

N = 7;

#debug_Sylvester_SMW_WEP(nx, nz, delta, N)




##################################################################################################################
nz = 11;
nx = nz + 4;
delta = 0.1;

println("\n     Test Tausch FD WEP-format\n")
nep_jar_wep = nep_gallery("waveguide", nx, nz, "Jarlebring", "fD", "weP", delta)

gmres_kwargs = ((:maxiter,187), (:restart,187), (:log,false))
function wep_gmres_linsolvercreator(nep::NEP, λ)
    return gmres_linsolvercreator(nep, λ, gmres_kwargs)
end


λ_jar_wep=NaN;
x_jar_wep=NaN
try
    λ_jar_wep,x_jar_wep =res_inv(nep_jar_wep, displaylevel=1, λ=-0.5-0.4im, maxit = 50, tolerance = 1e-10, v=ones(Complex128,nx*nz+2*nz), linsolvercreator=wep_gmres_linsolvercreator)
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ_jar_wep = e.λ
    x_jar_wep = e.v
end
println("Resnorm: ", compute_resnorm(nep_jar_wep, λ_jar_wep, x_jar_wep))
println("Eigenvalue: ", λ_jar_wep)
println("Eigenvector norm: ", norm(x_jar_wep), "\n")



##################################################################################################################
nz = 105;
nx = nz + 4;
delta = 0.1;

println("\n     Test Tausch FD SPFM\n")
nep_tausch = nep_gallery("waveguide", nx, nz, "TAUSCH", "fD", "SPMF", delta)

λ_tausch=NaN;
x_tausch=NaN
try
    λ_tausch,x_tausch =res_inv(nep_tausch;displaylevel=1, λ=-0.015-5.1im, maxit = 50, tolerance = 1e-10, v=ones(Complex128,nx*nz+2*nz))
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ_tausch=e.λ
    x_tausch=e.v
end
println("Resnorm: ", compute_resnorm(nep_tausch,λ_tausch,x_tausch))
println("Eigenvalue: ", λ_tausch)
println("Eigenvector norm: ", norm(x_tausch), "\n")


##################################################################################################################
nz = 105;
nx = nz + 4;
delta = 0.1;

println("\n     Test Jarlebring FD SPFM\n")
nep_jar = nep_gallery("waveguide", nx, nz, "jArleBRIng", "fD", "SpmF", delta)

λ_jar=NaN;
x_jar=NaN
try
    λ_jar,x_jar =res_inv(nep_jar;displaylevel=1, λ=-0.5-0.4im, maxit = 50, tolerance = 1e-10, v=ones(Complex128,nx*nz+2*nz))
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ_jar=e.λ
    x_jar=e.v
end
println("Resnorm: ", compute_resnorm(nep_jar,λ_jar,x_jar))
println("Eigenvalue: ", λ_jar)
println("Eigenvector norm: ", norm(x_jar), "\n")


##################################################################################################################
delta = 0.1;

println("\n     Different DEBUG-tests\n")

matlab_debug_WEP_FD(119, 115, delta)

matlab_debug_full_matrix_WEP_FD_SPMF(21, 17, delta)

fft_debug_mateq(431, 427, delta)

debug_sqrtm_schur(281)

debug_sqrt_derivative()

debug_Mlincomb_FD_WEP(29, 31, delta)
