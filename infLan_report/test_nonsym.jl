using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot, DelimitedFiles
import ..NEPSolver.ilan;
import ..NEPSolver.iar;

function fix_err_vec!(v)
    for j=1:length(v)-1
        if v[j+1]>v[j]
            v[j+1]=v[j]
        end
    end
end


include("../src/method_ilan.jl");
include("../src/method_iar.jl");


n=500
m=250
Random.seed!(1) # reset the random seed
# K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n];
# A1 = sparse(K, J, rand(3*n-2))+sparse(K, J, rand(3*n-2))*im;
# A2 = sparse(K, J, rand(3*n-2))+sparse(K, J, rand(3*n-2))*im;
# A3 = sparse(K, J, rand(3*n-2))+sparse(K, J, rand(3*n-2))*im;
# A4 = sparse(K, J, rand(3*n-2))+sparse(K, J, rand(3*n-2))*im;
A1=spdiagm(-1 => -ones(n-1), 0 => -2*ones(n), 1 => -ones(n-1))*n
A2=one(A1)
A3=spdiagm(-1 => ones(n-1), 0 => 0*ones(n) , 1 => ones(n-1))
A4=spdiagm(-1 =>  1im*ones(n-1) )




f1= S -> one(S)
f2= S -> -S
f3= S -> S*sin(S)
f4= S -> exp(-S)
nep1=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])

nA1=norm(A1); nA2=norm(A2); nA3=norm(A3); nA4=norm(A4);

rel_err=(λ,z)->compute_resnorm(nep1,λ,z)/((
nA1*abs(f1(λ))+
nA2*abs(f2(λ))+
nA3*abs(f3(λ))+
nA4*abs(f4(λ))
)*norm(z));
σ=0;
γ=1;
λ,_,err_iar=iar(nep1;σ=σ,γ=γ,Neig=100,displaylevel=1,maxit=m,tol=1e-12,check_error_every=1,errmeasure=rel_err)

O=zero(A1);
AA1=[O A1; transpose(A1) O];
AA2=[O A2; transpose(A2) O];
AA3=[O A3; transpose(A3) O];
AA4=[O A4; transpose(A4) O];

nep=SPMF_NEP([AA1,AA2,AA3,AA4],[f1,f2,f3,f4])
nep=DerSPMF(nep,σ,2*m)

V,H,ω,HH=ilan(nep;σ=σ,γ=γ,Neig=200,displaylevel=1,maxit=200,tol=1e-12,check_error_every=1)
#V,_,_=svd(V)
Q=V;

# Create a projected NEP
m=size(V,2)
pnep=create_proj_NEP(nep,m); # maxsize=mm
set_projectmatrices!(pnep,V,V);

err_lifted=(λ,z)->rel_err(λ,(Q*z)[n+1:end])
λ1,_,err=iar(pnep;σ=σ,γ=γ,Neig=100,displaylevel=1,maxit=200,tol=1e-12,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);

# sort error
for j=1:100
    err[:,j]=sort(err[:,j];rev=true);
    err_iar[:,j]=sort(err_iar[:,j];rev=true);
    #fix_err_vec!(view(err,:,j))
    #fix_err_vec!(view(err_iar,:,j))
end

# plot err hist
figure()
m=100
q=1
for j=1:m
    semilogy(1:q:m,err[1:q:m,j],color="black",linestyle="-");
    semilogy(1:q:m,err_iar[1:q:m,j],color="red",linestyle="--");
end

ylim(ymax=1)

figure()
plot(real(λ1),imag(λ1),marker="o",markerfacecolor=:none,c=:red,linestyle=:none) # Ritz values
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none) # Ritz values

# save the output

# save in a file the eigenvalues and Ritz pairs
writedlm("eig_rand_unsymm_ilan.csv",[real(λ1) imag(λ1)],",")
writedlm("eig_rand_unsymm_iar.csv",[real(λ) imag(λ)],",")



# save the error-matrix that will be loaded in tikz
err_print=NaN*ones(m,m+1)
err_print[1:m,1]=1:m
err_print[1:m,2:m+1]=err[1:m,2:m+1]
writedlm("err_rand_unsymm_ilan.csv",err_print,",")

err_print=NaN*ones(m,m+1)
err_print[1:m,1]=1:m
err_print[1:m,2:m+1]=err_iar[1:m,2:m+1]
writedlm("err_rand_unsymm_iar.csv",err_iar,",")
