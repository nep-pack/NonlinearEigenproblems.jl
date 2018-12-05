using NonlinearEigenproblems, Random, SparseArrays, Revise, DelimitedFiles
import ..NEPSolver.ilan;
import ..NEPSolver.iar;
import ..NEPSolver.tiar;


include("../src/method_ilan.jl");
include("../src/method_iar.jl");
include("../src/method_tiar.jl");


# load the Voss symmetric DEP
n=320; nep=nep_gallery("dep_symm_double",n)
nA1=opnorm(nep.A[1],Inf)
nA2=opnorm(nep.A[2],Inf)
err_rel=(λ,z)->compute_resnorm(nep,λ,z)/(abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[1]))*nA2);
_,_,_,_,conv_eig_hist_ilan=ilan(Float64,nep;Neig=200,displaylevel=1,maxit=100,tol=1e-12,check_error_every=10,errmeasure=err_rel)
_,_,_,_,conv_eig_hist_tiar=tiar(Float64,nep;Neig=200,displaylevel=1,maxit=100,tol=1e-12,check_error_every=10,errmeasure=err_rel)

for i=1:length(conv_eig_hist_ilan)
  println(conv_eig_hist_ilan[i])
end

for i=1:length(conv_eig_hist_tiar)
  println(conv_eig_hist_tiar[i])
end

# now export the conv-matrix that will be loaded in tikz
m=100
conv_eig_hist_tiar_print=ones(10,2)
conv_eig_hist_tiar_print[1:10,1]=10:10:100
conv_eig_hist_tiar_print[:,2]=conv_eig_hist_tiar[10:10:m]
writedlm("conv_eig_hist_tiar_print.csv",conv_eig_hist_tiar_print,",")

conv_eig_hist_ilan_print=ones(10,2)
conv_eig_hist_ilan_print[1:10,1]=10:10:100
conv_eig_hist_ilan_print[:,2]=conv_eig_hist_ilan[10:10:m]
writedlm("conv_eig_hist_ilan_print.csv",conv_eig_hist_ilan_print,",")
