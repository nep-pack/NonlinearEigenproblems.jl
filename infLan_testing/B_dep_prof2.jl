# TODO: compare also with the original function that does not perform the lr approximation
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


n=100000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';

nep=DEP([A1,A2,A3],[0,1,0.8])
fv=get_fv(nep); p=length(fv);    Av=get_Av(nep)

k=300; T=Float64

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


# precompute derivatives and FDH matrices (TODO: move in a preomputation function)
SS=diagm(0 => σ*ones(T,2k+2)) + diagm(-1 => γ*(1:2k+1))
fD=zeros(T,2*k+2,p)
for t=1:p fD[:,t]=fv[t](SS)[:,1] end
FDH=Vector{Array{T,2}}(undef,p)
for t=1:p FDH[t]=zeros(T,k+1,k+1)
    for i=1:k+1 for j=1:k+1
        FDH[t][i,j]=fD[i+j,t];
    end end
end
function Bmult!(k,Z,Qn,Av,FDH,G)
    # B-multiplication
    Z[:,:]=zero(Z);
    @inbounds for t=1:length(Av)
        Z[:,:] .+= Av[t]*Qn*(G.*view(FDH[t],1:k+1,1:k+1));
    end
    return Z
end

function Bmult_lr2!(k,Z,Qn,Av,G,vv)
    # B-multiplication
    Z[:,:]=zero(Z)
    ZZ=zero(Z)  # aux matrix for pre-allocation
    # low-rank factorization of G
    tolG=1e-12; U,S,V=svd(G[1:k+1,1:k+1]);q=sum(S.>tolG*ones(length(S)))
    U=view(U,:,1:q).*sqrt.(S[1:q]')
    V=view(V,:,1:q).*sqrt.(S[1:q]');
    Z[:,1]=Qn[:,1] # first matrix: fix for different \sigma
    @inbounds for t=2:length(Av)
        #Z[:,:] += Av[t]*((Qn*(U.*view(vv,:,t-1)))*(V.*view(vv,:,t-1))')
        mul!(ZZ,Av[t],Qn*(U.*view(vv,:,t-1))*(V.*view(vv,:,t-1))'); Z[:,:] += ZZ
        #Z[:,:] = muladd(Av[t],((Qn*(U.*view(vv,:,t-1)))*(V.*view(vv,:,t-1))'),Z[:,:])
    end
    Z[:,:] = -Z[:,:]
    return Z
end


Qn=rand(T,n,k+1)
Z=rand(T,n,k+1)

@btime begin Z3=Bmult!(k,copy(Z),Qn,Av,FDH,G) end
@btime begin Z2=Bmult_lr2!(k,copy(Z),Qn,Av,G,vv) end
#@btime begin Z1=Bmult_lr!(k,copy(Z),Qn,Av,G,vv) end


Z3=Bmult!(k,copy(Z),Qn,Av,FDH,G)
Z2=Bmult_lr2!(k,copy(Z),Qn,Av,G,vv)
#Z1=Bmult_lr!(k,copy(Z),Qn,Av,G,vv)

#display(norm(Z1-Z2))
#display(norm(Z1-Z3))
display(norm(Z2-Z3))
