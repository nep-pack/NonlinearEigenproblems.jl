#Methods to transform a PEP to a companion linearized form and solve the corresponding linearized pencil
using SparseArrays

export companion
export polyeig #Wrapper around the solver for a linearized PEP pencil


"""
    E,A = companion(nep::Pep);

Linearizes a  polynomial eigenvalue problem (PEP) a to the companion form, as in the paper by Mehrmann and Voss.
More precisely, for a k-th degree PEP with n-by-n coefficient matrices,
this returns matrices E and A, both kn-by-kn, corresponding to the linearized problem
```math
Ax = λEx
```

# Example
```julia-repl
julia> pep = nep_gallery("pep0");
julia> E,A = companion(pep);
julia> λ, V = eigen(A,E);
julia> minimum(svd(compute_Mder(pep,λ[1])).S)
2.703104679937224e-12
```

# References
* V. Mehrmann and H. Voss, Non-linear eigenvalue problems, a challenge for modern eigenvalue methods, GAMM‐Mitteilungen (2004)
"""
    function companion(pep::PEP)

        n = size(pep,1) #Size of monomial coefficient matrices
        d = size(pep.A,1)-1 #Degree of pep
        T = eltype(pep.A[1]) #Deduce the type of elements in the PEP matrices

        ##### Check sparsity of the problem and allocate memory accordingly #####
        if issparse(pep)
            E = spzeros(T, d*n, d*n)
            A = spzeros(T, n*d, n*d)

            #(d-1)n-by-(d-1)n matrix (Used to construct both E and A)
            Iblock = sparse(kron(Matrix{T}(I, d-1, d-1), Matrix{T}(I, n, n)))
        else
            E = zeros(T, d*n, d*n)
            A = zeros(T, n*d, n*d)

            Iblock = kron(Matrix{T}(I, d-1, d-1), Matrix{T}(I, n, n))
        end

        E[1:n,1:n] = pep.A[d+1];#Fill block (1,1)
        E[n+1:d*n,n+1:d*n] = Iblock;#Fill all blocks on the diagonal with eye(n)

        #####Construct A #####

        #First row block of A
        for i=1:d
           A[1:n,(i-1)*n+1:i*n] = pep.A[d-i+1];
        end
        #Lower part of A
        A[n+1:d*n,1:(d-1)*n] = T(-1.0)*Iblock

        return E,-A


    end


"""
    λ,v = polyeig([eltype],nep::PEP,[eigsolvertype,])

Linearizes a  polynomial eigenvalue problem (PEP) a to the companion form
and solves the corresponding linear eigenvalue problem; see [`companion`](@ref).
The `eigsolvertype` is optinal can be used to specify how the linear problem
is solved; see [`eig_solve`](@ref), and [`EigSolver`](@ref).

# Example
```julia-repl
julia> pep = nep_gallery("pep0");
julia> λ,V = polyeig(pep);
julia> minimum(svd(compute_Mder(pep,λ[1])).S)
2.1724582040065456e-14
julia> norm(compute_Mlincomb(pep,λ[2],vec(V[:,2])))
1.2210363164200074e-12
```
"""
    polyeig(pep::PEP,vargs...)=polyeig(ComplexF64,pep,vargs...)
    function polyeig(::Type{T},pep::PEP,eigsolvertype::Type=DefaultEigSolver) where T

        #Linearize to Ax = λEx
        E,A = companion(pep);

        solver::EigSolver = eigsolvertype(A,E);
        D,V = eig_solve(solver,target=one(T),nev=size(A,1));

        D = Vector{T}(D)
        V = Matrix{T}(V)

        return D,V[1:size(pep,1),:]
    end



"""
    λ,v = polyeig([eltype],nep::ChebPEP)

Computes a companion linearization for the NEP represented
in a Chebyshev basis, and returns eigenpairs.

# Example
```julia
julia> using LinearAlgebra
julia> nep=nep_gallery("dep0");
julia> chebpep=ChebPEP(nep,9,-3,1,cosine_formula_cutoff=5);
julia> (λv,V)=polyeig(chebpep);
julia> ii=argmin(abs.(λv));
julia> λ=λv[ii];
julia> v=V[:,ii];
julia> norm(compute_Mlincomb(chebpep,λ,v))
1.3543968603949142e-14
julia> # Actually, it's not a bad approx to the original NEP either
julia> norm(compute_Mlincomb(nep,λ,v))
4.326355966047557e-6
```

# References:

* Amiraslani, A., Corless, R. M. & Lancaster, P. "Linearization of matrix polynomials expressed in poly-nomial bases" IMA J. Numer. Anal.,29 (2009): 141–157.
* Effenberger and Kressner. "Chebyshev interpolation for nonlinear eigenvalue problems." BIT Numerical Mathematics 52.4 (2012): 933-951.

"""
polyeig(pep::ChebPEP,vargs...)=polyeig(ComplexF64,pep,vargs...)
function polyeig(::Type{T}, pep::ChebPEP{TT,Ftype}) where {T,Ftype,TT}
    k=pep.k
    Fk=get_Av(pep);

    n=size(pep,1);
    L0=zeros(T,n*(k-1),n*(k-1));
    L1=zeros(T,n*(k-1),n*(k-1));
    II=Matrix{Ftype}(I,n,n);

    for j=1:(k-2)
        L0[((j-1)*n) .+ (1:n), j*n .+ (1:n)]=II
        L0[j*n .+ (1:n), ((j-1)*n) .+ (1:n)]=II
    end

    for j=1:k-1
        L0[((k-2)*n) .+ (1:n), (n*(j-1)).+ (1:n)]=-Fk[j];
    end

    L0[((k-2)*n) .+ (1:n), (n*(k-3)) .+ (1:n)] += Fk[k];

    for j=1:k-2
        factor=2;
        if (j==1)
            factor=1
        end
        L1[((j-1)*n) .+ (1:n), ((j-1)*n) .+ (1:n)]=factor*II
    end
    L1[((k-2)*n) .+ (1:n),
       ((k-2)*n) .+ (1:n)]=2*Fk[k]

    LL=eigen(L0,L1);

    a=pep.a; b=pep.b;
    return ((b-a)*(LL.values .+ 1)/2 .+ a,
            hcat(map(i -> normalize(LL.vectors[1:n,i]), 1:(n*(k-1)))...))

end
