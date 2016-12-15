#  This is the first code in NEP-pack


println("hello world from NEP-pack")


n=5;

srand(0) # reset the random seed
A0=randn(n,n);
A1=randn(n,n);
I=eye(n,n);
tau=1;

# short-hand definition of functions
M(λ)=-λ*I+A0+A1*exp(-tau*λ)
Mp(λ)=-I-tau*A1*exp(-tau*λ)


v=randn(n,1);
c=v; λ=0;
for k=1:10
    println("Iteration:",k," resnorm:",norm(M(λ)*v))

    # Note: Julia has no comma in matrices
    J=[M(λ) Mp(λ)*v; c' 0];
    F=[M(λ)*v; c'*v-1];
    delta=-J\F;
    # Note: Submatrix indexing with [] not with ()
    v=v+delta[1:n];
    λ=λ+delta[n+1];
    
end


