using LinearAlgebra, Revise
function Gprod_test()

    n=5
    m=4
    # precompute the symmetrizer coefficients
    G=zeros(m+1,m+1);
    for i=1:m+1 G[i,1]=1/i end
    for j=1:m
        for i=1:m+1
            G[i,j+1]=(G[i,j]*j)/(i+j);
        end
    end
    tolG=1e-14
    U,S,V=svd(G);
    q=sum(S.>tolG*ones(length(S)))
    U=U[:,1:q]*diagm(0=>sqrt.(S[1:q]));
    V=V[:,1:q]*diagm(0=>sqrt.(S[1:q]));
    display(norm(G-U*V'))

    FDH=rand(m+1,m+1);
    Z=zeros(n,m+1);
    Q=rand(n,m+1);
    A=rand(n,n)
    Z=A*Q*(G.*FDH)

    ZZ=zeros(n,m+1)
    for j=1:q
        ZZ=ZZ+A*(broadcast(*,broadcast(*,Q,view(U,:,j)')*FDH,view(V,:,j)'));
    end
    display(norm(Z-ZZ))
end

Gprod_test()
