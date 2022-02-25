# A delay eigenvalue problem
function dep0(n::Int=5)
    msws_rng = MSWS_RNG()
    A0 = gen_rng_mat(msws_rng,n,n)
    A1 = gen_rng_mat(msws_rng,n,n)
    tau = 1.0
    nep = DEP([A0,A1],[0,tau])
    return nep
end


# A delay eigenvalue problem
function dep2(n::Int=5)
    msws_rng = MSWS_RNG()
    A0 = gen_rng_mat(msws_rng,n,n)
    A1 = gen_rng_mat(msws_rng,n,n)
    A2 = gen_rng_mat(msws_rng,n,n)
    tau = 1.1
    tau2 = 0.4
    nep = DEP([A0,A1,A2],[0,tau,tau2])
    return nep
end


# A delay eigenvalue problem with sparse matrices
function dep0_sparse(n::Integer=100, p::Real=0.25)
    msws_rng = MSWS_RNG()
    A0 = sparse(1:n,1:n,vec(gen_rng_mat(msws_rng,n,1)))+gen_rng_spmat(msws_rng,n,n,p)
    A1 = sparse(1:n,1:n,vec(gen_rng_mat(msws_rng,n,1)))+gen_rng_spmat(msws_rng,n,n,p)
    tau = 1.0
    nep = DEP([A0,A1],[0,tau])
    return nep
end

# A delay eigenvalue problem with sparse matrices
function dep0_sparse_symmetric(n::Integer=100, p::Real=0.15)
    msws_rng = MSWS_RNG()
    A0 = sparse(1:n,1:n,vec(gen_rng_mat(msws_rng,n,1)))+gen_rng_spmat(msws_rng,n,n,p)
    A0 = A0+A0'
    A1 = sparse(1:n,1:n,vec(gen_rng_mat(msws_rng,n,1)))+gen_rng_spmat(msws_rng,n,n,p)
    A1 = A1+A1'
    A2 = sparse(1:n,1:n,vec(gen_rng_mat(msws_rng,n,1)))+gen_rng_spmat(msws_rng,n,n,p)
    A2 = A2+A2';
    tau = 1.0
    tau2 = 0.8
    nep = DEP([A0,A1,A2],[0,tau,tau2])
    return nep
end


# A delay eigenvalue problem with sparse tridiagonal matrices
function dep0_tridiag(n::Integer=100)
    msws_rng = MSWS_RNG()
    K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n] # sparsity pattern of tridiag matrix
    A0 = sparse(K, J, vec(gen_rng_mat(msws_rng,3*n-2,1)))
    A1 = sparse(K, J, vec(gen_rng_mat(msws_rng,3*n-2,1)))
    tau = 1.0
    nep = DEP([A0,A1],[0,tau])
    return nep
end


# A polynomial eigenvalue problem
function pep0(n::Integer=200)
    msws_rng = MSWS_RNG()
    A0=gen_rng_mat(msws_rng,n,n)
    A1=gen_rng_mat(msws_rng,n,n)
    A2=gen_rng_mat(msws_rng,n,n)
    A=[A0,A1,A2]
    nep=PEP(A)
    return nep
end


# A polynomial eigenvalue problem with symmetric matrices
function pep0_sym(n::Integer=200)
    msws_rng = MSWS_RNG()
    A0 = Symmetric(gen_rng_mat(msws_rng,n,n))
    A1 = Symmetric(gen_rng_mat(msws_rng,n,n))
    A2 = Symmetric(gen_rng_mat(msws_rng,n,n))
    A = [A0, A1, A2]
    nep = PEP(A)
    return nep
end


# A sparse polynomial eigenvalue problem
function pep0_sparse(n::Integer=200, p::Real=0.03)
    msws_rng = MSWS_RNG()
    A0=gen_rng_spmat(msws_rng,n,n,p)
    A1=gen_rng_spmat(msws_rng,n,n,p)
    A2=gen_rng_spmat(msws_rng,n,n,p)
    A=[A0,A1,A2]
    nep=PEP(A)
    return nep
end


# Helper for random matrices, stable over releases
# B. Widynski, Middle Square Weyl Sequence RNG, arXiv 1704.00358
mutable struct MSWS_RNG
    x::UInt128
    w::UInt128
    s::UInt128
    function MSWS_RNG(seed::UInt128 = zero(UInt128))
        base = 0x9ef09a97ac0f9ecaef01c4f2db0958c9
        s = (seed << 1) + base
        x = 0x1de568e1a1ca1b593cbf13f7407cf43e
        w = 0xd4ac5c288559e14a5fafc1b7df9f9e0e
        return new(x, w, s)
    end
end

function gen_rng_int(rng::MSWS_RNG)
    rng.x *= rng.x
    rng.x += (rng.w += rng.s)
    rng.x = (rng.x>>64) | (rng.x<<64)
    return UInt64(rng.x & typemax(UInt64))
end

function gen_rng_float(rng::MSWS_RNG)
    return Float64(gen_rng_int(rng)/typemax(UInt64))
end

function gen_rng_mat(rng::MSWS_RNG, n::Integer, m::Integer)
    A = zeros(Float64,n,m)
    for c = 1:m
        for r = 1:n
            A[r,c] = 1-2*gen_rng_float(rng)
        end
    end
    return A
end

function gen_rng_spmat(rng::MSWS_RNG, n::Integer, m::Integer, p::Real)
    n=Int64(n)
    m=Int64(m)
    nonzeros = round(p*m*n)
    dict = Dict{Tuple{Int64,Int64},Float64}()
    for i = 1:nonzeros
        r = mod(gen_rng_int(rng),n) + 1
        c = mod(gen_rng_int(rng),m) + 1
        dict[r,c] = 1-2*gen_rng_float(rng)
    end
    idxes = collect(keys(dict))
    nnzs = length(idxes)
    r = zeros(Int64,nnzs)
    c = zeros(Int64,nnzs)
    vals = zeros(Float64,nnzs)
    for i = 1:nnzs
        r[i] = idxes[i][1]
        c[i] = idxes[i][2]
        vals[i] = dict[idxes[i]]
    end
    return sparse(r,c,vals,n,m)
end
