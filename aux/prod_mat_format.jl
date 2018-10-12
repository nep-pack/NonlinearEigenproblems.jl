using Random

n=2000;
k=100;
V=rand(n,k)+rand(n,k)*im;
W=rand(n,k)+rand(n,k)*im;
α=V[:]⋅W[:]
# this is equivalent to V[:]⋅W[:]
β=sum(sum(conj(V).*W,dims=1));
norm(α-β)
