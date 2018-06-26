n=10
L=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)
L=n^2*kron(L,L)

x = linspace(0,pi,n)
b=broadcast((x,y)->100*abs(sin(x+y)),x,x.')
a=broadcast((x,y)->-sin(x)*sin(y),x,x.')
B=sparse(1:n^2,1:n^2,b[:])
A=L+sparse(1:n^2,1:n^2,a[:])
