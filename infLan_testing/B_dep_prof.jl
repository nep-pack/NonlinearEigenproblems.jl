using NonlinearEigenproblems, Random, SparseArrays, Revise, LinearAlgebra, BenchmarkTools



function Bmult_lr!(k,Z,Qn,Av,G,vv)
    # B-multiplication
    Z[:,:]=zero(Z);
    # low-rank factorization of G
    tolG=1e-12; Ug,Sg,Vg=svd(G[1:k+1,1:k+1]);q=sum(Sg.>tolG*ones(length(Sg)))
    Ug=broadcast(*,view(Ug,:,1:q),sqrt.(Sg[1:q])');
    Vg=broadcast(*,view(Vg,:,1:q),sqrt.(Sg[1:q])');
    Z[:,1]=-Qn[:,1] # first matrix: fix for different \sigma
    # TODO: avoid materialization
    @inbounds for t=2:length(Av)
        #@simd
        for j=1:q
            ZZ=Qn.*(Ug[:,j]')
            ZZ=-broadcast(*,ZZ*vv[:,t-1],vv[:,t-1]')
            Z[:,:] .+= Av[t]*(ZZ.*(Vg[:,j]'));
        end
    end
    return Z
end


n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';

nep=DEP([A1,A2,A3],[0,1,0.8])
fv=get_fv(nep); p=length(fv);    Av=get_Av(nep)

k=100; T=Float64

# precomputation for DEP (TODO: move in a preomputation function)
vv=zeros(T,k+1,p-1)
τ=nep.tauv
γ=1; σ=0
for j=1:p-1
    vv[:,j]=sqrt(τ[j]*γ)*exp(-σ)*(-τ[j]*γ).^(0:k)
end

# precompute the symmetrizer coefficients
G=zeros(T,k+1,k+1);
for i=1:k+1 G[i,1]=1/i end
for j=1:k
    for i=1:k+1
        G[i,j+1]=(G[i,j]*j)/(i+j);
    end
end

Qn=rand(T,n,k+1)
Z=rand(T,n,k+1)
@btime Z1=Bmult_lr!(k,copy(Z),Qn,Av,G,vv)



function Bmult_lr2!(k,Z,Qn,Av,G,vv)
    # B-multiplication
    Z[:,:]=zero(Z);
    # low-rank factorization of G
    tolG=1e-12; U,S,V=svd(G[1:k+1,1:k+1]);q=sum(S.>tolG*ones(length(S)))
    U=broadcast(*,view(U,:,1:q),sqrt.(S[1:q])');
    V=broadcast(*,view(V,:,1:q),sqrt.(S[1:q])');
    Z[:,1]=-Qn[:,1] # first matrix: fix for different \sigma
    @inbounds for t=2:length(Av)
        Z[:,:] -= Av[t]*((Qn*(U.*view(vv,:,t-1)))*(V.*view(vv,:,t-1))')
    end
    return Z
end


@btime Z2=Bmult_lr2!(k,copy(Z),Qn,Av,G,vv)

display(norm(Z1-Z2))
