#  This is the first code in NEP-pack


println("hello world from NEP-pack")


n=5;


A0=randn(n,n);
A1=randn(n,n);
I=eye(n,n);
tau=1;
# short-hand definition of functions
M(lambda)=-lambda*I+A0+A1*exp(-tau*lambda)
Mp(lambda)=-I-tau*A1*exp(-tau*lambda)


v=randn(n,1);
c=v; lambda=0;
for k=1:10
    println("Iteration:",k," resnorm:",norm(M(lambda)*v))

    # Note: Julia has no comma in matrices
    J=[M(lambda) Mp(lambda)*v; c' 0];
    F=[M(lambda)*v; c'*v-1];
    delta=-J\F;
    # Note: Submatrix indexing with [] not with ()
    v=v+delta[1:n];
    lambda=lambda+delta[n+1];
    
end


