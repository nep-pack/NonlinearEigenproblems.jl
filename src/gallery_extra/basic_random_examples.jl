# A delay eigenvalue problem
function dep0(n::Int=5)
    Random.seed!(0) # reset the random seed
    A0 = randn(n,n)
    A1 = randn(n,n)
    tau = 1.0
    nep = DEP([A0,A1],[0,tau])
    return nep
end


# A delay eigenvalue problem with sparse matrices
function dep0_sparse(n::Integer=100, p::Real=0.25)
    Random.seed!(0) # reset the random seed
    A0 = sparse(1:n,1:n,rand(n))+sprand(n,n,p)
    A1 = sparse(1:n,1:n,rand(n))+sprand(n,n,p)
    tau = 1.0
    nep = DEP([A0,A1],[0,tau])
    return nep
end


# A delay eigenvalue problem with sparse tridiagonal matrices
function dep0_tridiag(n::Integer=100)
    Random.seed!(1) # reset the random seed
    K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n] # sparsity pattern of tridiag matrix
    A0 = sparse(K, J, rand(3*n-2))
    A1 = sparse(K, J, rand(3*n-2))
    tau = 1.0
    nep = DEP([A0,A1],[0,tau])
    return nep
end


# A polynomial eigenvalue problem
function pep0(n::Integer=200)
    Random.seed!(0)
    A0=randn(n,n)
    A1=randn(n,n)
    A2=randn(n,n)
    A=[A0,A1,A2]
    nep=PEP(A)
    return nep
end


# A polynomial eigenvalue problem with symmetric matrices
function pep0_sym(n::Integer=200)
    Random.seed!(0)
    A0 = Symmetric(randn(n,n))
    A1 = Symmetric(randn(n,n))
    A2 = Symmetric(randn(n,n))
    A = [A0]#, A1, A2]
    nep = PEP(A)
    return nep
end


# A sparse polynomial eigenvalue problem
function pep0_sparse(n::Integer=200, p::Real=0.03)
    Random.seed!(0)
    A0=sprandn(n,n,p)
    A1=sprandn(n,n,p)
    A2=sprandn(n,n,p)
    A=[A0,A1,A2]
    nep=PEP(A)
    return nep
end
