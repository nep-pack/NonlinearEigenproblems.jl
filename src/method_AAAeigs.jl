import SparseArrays: findnz


# Struct holding information for invoking svAAA in AAAPencil
struct AAACorkLinearization{T<:Number} <: CorkLinearization
    Z::AbstractArray{T}
    mmax::Int
    tol::Real
    cleanup::Bool
    tol_cln::Real
    return_details::Bool
    logger::Logger
    weighted::Bool
    function AAACorkLinearization(
            Z::AbstractArray{T};
            mmax::Int=100,
            tol::Real=eps(real(T)) * 1e3,
            cleanup::Bool=true,
            tol_cln::Real = min(eps(real(T)),tol),
            return_details::Bool=false,
            logger=0,
            weighted::Bool=false) where {T<:Number}

        @parse_logger_param!(logger)
        return new{T}(Z, mmax, tol, cleanup, tol_cln, return_details, logger, weighted)
    end
end

# Struct representing compact version of the AAA-pencil
struct AAAPencil{S<:AbstractVector{<:AbstractMatrix{<:Number}}}
    d::Int # Degree of polynomial part
    s::Int # Number of non-polynomial matrices and functions
    m::Int # Degree of AAA rational approximation
    PPCC::S # Vector containing all non-zero coefficient matrices of the nep
    ppff::AbstractVector{Function} # Vector containing all corresponding scalar functions of the nep
    NNZ::AbstractVector{Int} # Vector containing the degrees of the polynomial terms in 'PPCC'/'ppff'
    compactA::AbstractMatrix # Compact representation of L
    compactB::AbstractMatrix # Compact representation of B
end

# Constructor for AAAPencil using the svAAA-function
function AAAPencil(nep::NEP, is::AAACorkLinearization{T}) where T<:Number
    # If no access to coefficient matrices and scalar functions, AAA-approximation is not possible
    if !isa(nep, AbstractSPMF)
        error("AAAPencil: The given nep must be a sum of products of matrices and functions (i.e. `AbstractSPMF`)")
    end

    # Polynomial eigenvalue problem
    if isa(nep, PEP)
        Av = get_Av(nep)
        d = length(Av)-1
        if d < 2 # Linear eigenvalue problem
            error("AAAPencil: The given nep is a simple Av=λBv, solve with an appropriate linear eigensolver (e.g. `eigen(L,B)`)")
        end
        NNZ = findall(x->!x, map(iszero,Av))
        while NNZ[end] != d+1
            pop!(NNZ)
            d = d-1
        end
        PPCC = Av[NNZ]
        ppff = get_fv(nep)[NNZ]
        NNZ .-= 1
        s = 0
        m = 0
        z = Vector{T}(); fz = Matrix{T}(undef,0,0); ω = Vector{T}()
        err = Vector{real(T)}(); pol = Vector{T}(); rsd = Vector{T}(); zer = Vector{T}()

    # Polynomial + nonlinear eigenvalue problem
    elseif isa(nep, SPMFSumNEP{PEP,S} where S<:AbstractSPMF)
        Av_p = get_Av(nep.nep1)
        Av_f = get_Av(nep.nep2)
        d = length(Av_p)-1
        NNZ = findall(x->!x, map(iszero,Av_p))
        while NNZ[end] != d+1
            pop!(NNZ)
            d = d-1
        end
        PPCC = [Av_p[NNZ]             ; Av_f]
        ppff = [get_fv(nep.nep1)[NNZ] ; get_fv(nep.nep2)]
        NNZ .-= 1
        s = length(Av_f)
        z, fz, ω, err, pol, rsd, zer = svAAA(nep.nep2, is.Z,
            mmax=is.mmax, tol=is.tol, cleanup=is.cleanup, tol_cln=is.tol_cln, return_details=is.return_details, logger=is.logger, weighted=is.weighted)
        m = length(z)

    # Nonlinear + polynomial eigenvalue problem
    elseif isa(nep, SPMFSumNEP{S,PEP} where S<:AbstractSPMF)
        Av_p = get_Av(nep.nep2)
        Av_f = get_Av(nep.nep1)
        d = length(Av_p)-1
        NNZ = findall(x->!x, map(iszero,Av_p))
        while NNZ[end] != d+1
            pop!(NNZ)
            d = d-1
        end
        PPCC = [Av_p[NNZ]             ; Av_f]
        ppff = [get_fv(nep.nep2)[NNZ] ; get_fv(nep.nep1)]
        NNZ .-= 1
        s = length(Av_f)
        z, fz, ω, err, pol, rsd, zer = svAAA(nep.nep1, is.Z,
            mmax=is.mmax, tol=is.tol, cleanup=is.cleanup, tol_cln=is.tol_cln, return_details=is.return_details, logger=is.logger, weighted=is.weighted)
        m = length(z)

    # Completely nonlinear eigenvalue problem
    else
        NNZ = Vector{Int}()
        PPCC = get_Av(nep)
        ppff = get_fv(nep)
        d = 0
        s = length(PPCC)
        z, fz, ω, err, pol, rsd, zer = svAAA(nep, is.Z,
            mmax=is.mmax, tol=is.tol, cleanup=is.cleanup, tol_cln=is.tol_cln, return_details=is.return_details, logger=is.logger, weighted=is.weighted)
        m = length(z)
    end

    zfω = [z, fz, ω]
    prz = [pol, rsd, zer]

    # Form [P_A^T M^T] and [P_B^T N^T] using AAA-information
    compactA, compactB = get_compact_pencil(d, s, m, z, fz, ω, NNZ)
    push_info!(is.logger, "AAAPencil: Pencil is built with d=$d, s=$s and m=$m.")

    return AAAPencil(d, s, m, PPCC, ppff, NNZ, compactA, compactB), zfω, err, prz
end

# Get M and N from the linear relations and projections P_A and P_B
function get_compact_pencil(d::Int, s::Int, m::Int, z::AbstractVector{T}, fz::AbstractMatrix{T}, ω::AbstractVector{T}, NNZ::AbstractVector{Int}) where T<:Number
    # No polynomial part
    dt = length(NNZ)
    if dt == 0
        compactA = [ fz spdiagm(m,m-1, 0=>-ω[2:end].*z[1:end-1], -1=>ω[1:end-1].*z[2:end]) ]

        compactB = [ spzeros(m,s)  spdiagm(m,m-1, 0=>-ω[2:end], -1=>ω[1:end-1]) ]

    # Only constant term in polynomial part
    elseif d == 0
        compactA = [ spzeros(1,1+s+m) ;
                     spzeros(m,1) fz spdiagm(m,m-1, 0=> -ω[2:end].*z[1:end-1], -1=>ω[1:end-1].*z[2:end]) ones(m,1)]
        compactA[1,1] = 1 ; compactA[1,end] = -1

        compactB = [ spzeros(1,1+s+m) ;
                     spzeros(m,1+s)  spdiagm(m,m-1, 0=>-ω[2:end], -1=>ω[1:end-1]) spzeros(m,1)]

    # No nonlinear part
    elseif m == 0
        compactA = [ sparse(NNZ[1:end-1].+1,1:dt-1,ones(dt-1),d,dt-1) spzeros(d,1) spdiagm(d,d-1, -1=>ones(d-1))]

        compactB = [ spzeros(d,dt) spdiagm(d,d-1, 0=>ones(d-1))]
        compactB[d,dt] = -1

    # General case: polynomial + nonlinear part
    else
        compactA = [ sparse(NNZ[1:end-1].+1,1:dt-1,ones(dt-1),d,dt-1) spzeros(d,s+1) spdiagm(d,d-1, -1=>ones(d-1)) spzeros(d,m) ;
                     spzeros(m,dt)                                     fz            spzeros(m,d-1)                spdiagm(m,m-1, 0=> -ω[2:end].*z[1:end-1], -1=>ω[1:end-1].*z[2:end]) ones(m,1)]
        compactA[1,end] = -1 ;

        compactB = [ spzeros(d,dt+s) spdiagm(d,d-1, 0=>ones(d-1)) spzeros(d,m) ;
                     spzeros(m,dt+s+d-1)                          spdiagm(m,m-1, 0=>-ω[2:end], -1=>ω[1:end-1]) spzeros(m,1)]
        compactB[d,dt] = -1
    end

    return compactA, compactB
end


"""
    AAAeigs([eltype,] nep::NEP, Z::AbstractArray{T}; [logger=0,][mmax=100,][neigs=6,][maxit=min(max(10*neigs,30),100),][shifts=Vector{T}(),][linsolvercreator=FactorizeLinSolverCreator(max_factorizations=min(length(unique(shifts)),10)),][tol=eps(real(T))*1e6,][tol_appr=eps(real(T))*1e3,][v0=Vector{T}(),][errmeasure=ResidualErrmeasure(nep),][weighted=false,][cleanup_appr=true,][tol_cln=min(eps(real(T)),tol_appr),][return_details=false,][check_error_every=10,][inner_logger=0])

Find a few eigenvalues and eigenvectors of a nonlinear eigenvalue problem, using the `AAAeigs` algorithm.

# Arguments
- `nep`: An instance of a nonlinear eigenvalue problem.
- `Z`: An array containing the sample domain of the nep, given as a series of points.
- `mmax`: Max degree of AAA-approximation.
- `neigs`: Number of eigenvalues to find.
- `maxit`: Max number of total iterations.
- `shifts`: L vector of shifts used during the Krylov routine (empty -> shift=0)
- `tol`: Tolerance for residual.
- `tol_appr`: Tolerance for convergence of AAA-approximation.
- `v0`: Starting vector.
- `weighted`: Wether to use Weighted or Set-Valued AAA.
- `cleanup_appr`: Whether to detect spurious poles in AAA.
- `tol_cln`: Tolerance for cleanup of spurious poles.
- `return_details`: Whether to return solution details (see `AAASolutionDetails`).
- `check_error_every`: Check for convergence / termination every this number of iterations.
- `inner_logger`: Works in the same way as logger but for the AAA rational approximation (see `svAAA`)

See [`augnewton`](@ref) for other parameters.


# Return values
- `Λ`: Vector of eigenvalues of the nonlinear eigenvalue problem inside the sample set Z.
- `X`: Corresponding matrix of eigenvectors.
- `res`: Corresponding residuals.
- `details`: Solution details, if requested (see AAASolutionDetails).

# Example
```julia-repl
julia> Random.seed!(2022);
julia> nep=nep_gallery("nlevp_native_gun");
julia> m=250^2; r=300^2-200^2;
julia> Z=m-r.+2*r*rand(1000)+im*(rand(1000)*r.+sqrt(1e-13));
julia> Z=Z[(~).(imag.(Z).>sin.(acos.((real.(Z).-m)/r))*r.-sqrt(1e-13))];
julia> Z=[transpose([LinRange(m-r+1e-2,m+r-1e-2,250)' m-r.+2*r*(exp.(im*LinRange(0,pi,250)')/2 .+.5)]); Z[1:500]];
julia> shifts=r*[2/3,(1+im)/3,0,(-1+im)/3,-2/3].+m;
julia> λ,X=AAAeigs(nep,Z,shifts=shifts);
julia> [norm(compute_Mlincomb(nep,λ[i],X[:,i]))/norm(X[:,i]) for i=1:length(λ)]
6-element Vector{Float64}:
 5.282105614825289e-12
 1.675900913610527e-11
 4.971723316733914e-11
 9.437128735632508e-11
 1.2277986011104446e-10
 1.4546362158344052e-10
```
# References
- P. Lietaert, J. Pérez, B. Vandereycken, and K. Meerbergen. Automatic rational
  approximation and linearization of nonlinear eigenvalue problems. IMA journal
  of numerical analysis, 2018.
- R. Van Beeumen, K. Meerbergen, and W. Michiels. Compact rational Krylov
  methods for nonlinear eigenvalue problems. SIAM Journal on Matrix Analysis
  and Applications, 36(2):820-838, 2015.
"""
AAAeigs(nep::NEP, Z::AbstractArray{T}; params...) where T<:Number =
    AAAeigs(T, nep, Z; params...)
function AAAeigs(
        ::Type{T},
        nep::NEP,
        Z::AbstractArray{T};
        logger = 0,
        mmax::Int = 100,
        neigs::Real = 6,
        maxit::Int = min(max(10*neigs,30),100),
        shifts::Vector{T} = Vector{T}(),
        linsolvercreator = FactorizeLinSolverCreator(max_factorizations=min(length(unique(shifts)),10)),
        tol::Real = eps(real(T))*1e6,
        tol_appr::Real = eps(real(T))*1e3,
        v0::Vector{T} = Vector{T}(),
        errmeasure::ErrmeasureType = ResidualErrmeasure(nep),
        weighted::Bool = false,
        cleanup_appr::Bool = true,
        tol_cln::Real=min(eps(real(T)),tol_appr),
        return_details::Bool = false,
        check_error_every::Int = 10,
        inner_logger = 0) where T<:Number

    # Work out the inputs
    @parse_logger_param!(logger)
    @parse_logger_param!(inner_logger)

    n = size(nep,1)
    rmax = maxit
    jmax = maxit
    # Handle the shifts
    if isempty(shifts)
        shifts = [0.]
    end
    σ = Vector{T}(undef,maxit)
    mxtdivs = fld(maxit,length(shifts))
    mxtmods = mod(maxit,length(shifts))
    σ[1:length(shifts)*mxtdivs] = repeat(shifts,mxtdivs)
    σ[length(shifts)*mxtdivs+1:end] = shifts[1:mxtmods]

    # Perform AAA and obtain pencil
    is = AAACorkLinearization(Z, mmax=mmax, tol=tol_appr, weighted=weighted, cleanup=cleanup_appr, tol_cln=tol_cln, return_details=return_details, logger=inner_logger)
    L, zfω, err_appr, prz = AAAPencil(nep, is)
    #   Set (combinations of) size parameters
    d = L.d
    dt = length(L.NNZ)
    m = L.m
    k = d + m
    if d == 0 && dt != 0;
        k = k + 1; # d = 0/1 involves the same size linearization
    end
    s = L.s
    l = dt + s

    # Initializations
    #   Extended '(M-λN)'-matrix and dictionary for factorizations
    max_factorizations = min(length(unique(shifts)),10)
    MλN_factorizations = Dict{T, Any}()
    MλNext = spzeros(T,k,k)
    #   Intermediate matrices/vectors
    Y = zeros(T,k,l)
    Uhat = 0 # Otherwise scope Uhat in try-catch is wrong
    qq = Vector{T}(undef, n)
    #   Initialize 'growing' matrices of CORK quadruple
    if length(v0) != n # Case when not supplied or when wrong length
        v0 = randn(T,n)
    end
    Q = Matrix{T}(undef,n,min(50,rmax+1))
    Q[:,1] .= v0./norm(v0)
    U = zeros(T,rmax+1,k,jmax+1)
    u0 = [1; zeros(k-1)]
    U[1,:,1] .= u0./norm(u0)
    H = UpperHessenberg(Matrix{T}(undef,jmax+1,jmax))
    K = UpperHessenberg(Matrix{T}(undef,jmax+1,jmax))
    #   Matrices for eigenvalues and residues of all iterations
    if return_details
        Lam = Matrix{T}(undef,maxit,maxit)
        Res = Matrix{real(T)}(undef,maxit,maxit)
    end

    # Perform CORK iterations
    r = 1
    j = 1
    it = 1
    nconv = 0
    while (it<=maxit) && (nconv<neigs)
        # Level 1 orthogonalization
        #   Get factorization of MλNext
        if (σ[it] in keys(MλN_factorizations))
            MλNfact=MλN_factorizations[σ[it]]
        else
            MλNext = [ sparse(I,k,1) @views(L.compactA[:,l+1:end]-σ[it]*L.compactB[:,l+1:end])]
            MλNfact = factorize(MλNext)
            if  (length(keys(MλN_factorizations)) < max_factorizations )
                MλN_factorizations[σ[it]] = MλNfact
            end
        end
        #   Compute Y and solve system with MλNext
        Y .= @views ( σ[it].*L.compactB[:,1:l] .- L.compactA[:,1:l])
        Y = MλNfact\Y
        #   Compute new vector for Q using Uc and 'in-place'-summation
        u_c = @views ( U[1:r,1:k,j]*(L.compactB*[ sparse(I,l,l) ;
                                                  Y[2:end,:]    ]) )
        v1_hat = zeros(T, n)
        for i = 1:l
            mul!(qq, ( @view Q[:,1:r] ), ( @view u_c[:,i]) )
            mul!(v1_hat, L.PPCC[i], qq, one(T), one(T))
        end
        #   Solve system with NEP's matrix function
        solver = create_linsolver(linsolvercreator,nep,σ[it])
        v1_hat = lin_solve(solver,v1_hat)
        #   Gram-Schmidt-orthogonalization
        nv = norm(v1_hat)
        u1_hat = ( @view Q[:,1:r] )'*v1_hat
        mul!(v1_hat, ( @view Q[:,1:r] ), u1_hat, -one(T), one(T))
        #   DGKS-reorthogonalization
        ii = 0
        while (ii < 3) && (norm(v1_hat) < 1/sqrt(2)*nv)
            nv = norm(v1_hat)
            u1_new = ( @view Q[:,1:r] )'*v1_hat
            mul!(v1_hat, ( @view Q[:,1:r] ), u1_new, -one(T), one(T))
            u1_hat .+= u1_new
            ii = ii+1
        end
        #   expansion of Q-matrix
        nv = norm(v1_hat)
        if nv > eps()
            rnew = r + 1
            if size(Q,2) < rnew
                resized = zeros(T,n,min(2*r,rmax+1))
                resized[:,1:r] = Q
                Q = resized
            end
            Q[:,rnew] .= v1_hat./nv
            U[rnew,1:k,1:j] .= 0
            append!(u1_hat,nv)
        else
            rnew = r
        end

        # Level 2 orthogonalization
        #   Compute W and solve system with MλNext
        W = zeros(T,rnew,k)
        W .= u1_hat
        mul!( (@view W[:,2:end]), (@view U[1:rnew,1:k,j]), (@view L.compactB[:,l+1:end]) )
        try
            Uhat = W/MλNfact
        catch
            Uhat = transpose(transpose(MλNfact)\transpose(W))
        end
        #   Gram-Schmidt-orthogonalization
        U_rs = reshape(( @view U[1:rnew,:,1:j] ), rnew*k, j)
        uhat_rs = reshape(Uhat,rnew*k)
        nu = norm(uhat_rs)
        mul!((@view H[1:j,j]), U_rs', uhat_rs)
        mul!(uhat_rs, U_rs, (@view H[1:j,j]), -one(T), one(T))
        H[j+1,j] = norm(uhat_rs)
        #   DGKS-reorthogonalization
        ii = 0
        while (ii < 3) && (real(H[j+1,j]) < 1/sqrt(2)*nu)
            h_new = (U_rs')*uhat_rs
            mul!(uhat_rs, U_rs, h_new, -one(T), one(T))
            H[1:j,j] .= ( @view H[1:j,j] )  .+ h_new
            nu = real(H[j+1,j])
            H[j+1,j] = norm(uhat_rs)
            ii = ii+1
        end
        #   Update U, K and H from CORK quadruple
        U[1:rnew,:,j+1] .= Uhat./H[j+1,j]
        K[1:j,j] .= σ[it].*@view( H[1:j,j] )
        K[j,j] += 1 # Continuation vector contribution
        K[j+1,j] = H[j+1,j]*σ[it]

        # Check for convergence
        if(return_details || (rem(it,check_error_every)==0) || (it==maxit))
            # Compute Ritz pairs
            Λ, S = eigen!(K[1:j,1:j],H[1:j,1:j])
            X = @views( Q[:,1:rnew]*( U[1:rnew,1,1:j+1]*( H[1:j+1,1:j]*S ) ) )

            # Compute residues
            res = map(i -> estimate_error(errmeasure, Λ[i], (@view X[:,i])), 1:length(Λ))
            conv = abs.(res) .< tol
            nconv = isempty(conv) ? 0 : sum(conv)
            push_info!(logger, "AAAeigs iteration $it: $nconv of $it < $tol")
            push_iteration_info!(logger, 2, it, err=res, λ=Λ)
            idx = sortperm(res)

            # Save results in matrices
            if return_details
                Lam[:,it] .= T(NaN)       ; Lam[1:it,it] .= Λ[idx]
                Res[:,it] .= real(T)(NaN) ; Res[1:it,it] .= res[idx]
            end

            # Extract the converged Ritzpairs
            if (it==maxit) || (nconv >= neigs)
                nb_eigs = Int(min(length(Λ),neigs))
                Λ = Λ[idx[1:nb_eigs]]
                X = X[:,idx[1:nb_eigs]]
                res = res[idx[1:nb_eigs]]
            end
        end

        r = rnew
        j = j+1
        it = it+1
    end
    it = it-1

    # End of iterations
    if nconv < neigs && neigs != Inf
        idx = sortperm(res); res = res[idx]
        kk=length(Λ)
        Λ=Λ[idx[1:kk]]; X=X[:,idx[1:kk]]; res=res[1:kk]
        msg = "AAAeigs: Number of iterations exceeded. maxit=$maxit."
        throw(NoConvergenceException(Λ,Q,res,msg))
    end

    # Compile additional information
    details = AAASolutionDetails{T}()
    if return_details
        Lam = Lam[1:it,1:it]
        Res = Res[1:it,1:it]
        details = AAASolutionDetails(m, zfω, prz, err_appr, Lam, Res, it)
    end

    return Λ, X, res, details
end

# Sort given Ritz values Lam, according to their residue Res
function sort_eig!(Lam, Res)
    it = size(Res,2)
    for col = 1:it-1
        # Find the indices that sort the residues for the current column
        idx_sor = sortperm(@view Res[1:col,col])
        idx = Vector{Int}(1:col+1)
        idx_comb = zeros(Int,col+1)
        for row in idx_sor
            # Find the new ritz value that is closest to the previous in each row in order of smallest residue
            loc = argmin(abs.( @views (Lam[row,col].-Lam[idx,col+1] )))
            # Couple both entries in idx_comb
            idx_comb[row] = idx[loc]
            # The new ritz value that is selected is removed as index
            deleteat!(idx, findall(x->x.==idx[loc],idx))
        end
        # remaining new ritz value is placed at the end (not corresponding to any previous)
        idx_comb[end] = idx[end]
        # Sort the new column (of the ritz values and their residues)
        Lam[1:col+1,col+1] = Lam[idx_comb,col+1]
        Res[1:col+1,col+1] = Res[idx_comb,col+1]
    end
end

"""
    svAAA([eltype,] nep::AbstractSPMF, Z::AbstractArray{T}; [mmax=100,][tol=eps(real(T))*1e3,][cleanup=true,][tol_cln=min(eps(real(T),tol)),][return_details=false,][logger=0,][weighted=false])

Compute a Set-Valued AAA rational approximation of the (nonlinear) functions in the nonlinear eigenvalue problem.

# Arguments
- `nep`: An instance of a nonlinear eigenvalue problem.
- `Z`: An array containing the sample domain of the nep, given as a series of points.
- `mmax`: Max. number of support points/weights.
- `tol`: Relative tolerance for convergence.
- `cleanup`: Whether to detect spurious poles.
- `tol_cln`: Tolerance for cleanup of spurious poles.
- `return_details`: Whether to return details of the rational approximation.
- `logger`: Either a `Logger` object or an `Int`. If ii is an `Int`, a `PrintLogger(logger)` will be created with the appropriate display level.
- `weighted`: If true, the functions are weighted proportional to their appearance in the NEP; if false, the functions are scaled to be equal.

# Return values
- `z`: Vector of support points.
- `fz`: Matrix of function values in support points `z`.
- `ω`: Vector of weights of the barycentric representation.
- `err`: Vector of errors in successive iterations steps of svAAA.
- `pol`: Vector of poles of the rational approximation, if details requested.
- `rsd`: Vector of residues of the rational approximation, if details requested.
- `zer`: Vector of zeros of the rational approximation, if details requested.

# Example
```julia-repl
julia> Random.seed!(2022);
julia> nep=nep_gallery("nlevp_native_gun");
julia> typeof(nep)
SPMFSumNEP{PEP, SPMF_NEP{SparseMatrixCSC{Float64, Int64}, ComplexF64}}
julia> m=250^2; r=300^2-200^2;
julia> Z=m-r.+2*r*rand(1000)+im*(rand(1000)*r.+sqrt(1e-13));
julia> Z=Z[(~).(imag.(Z).>sin.(acos.((real.(Z).-m)/r))*r.-sqrt(1e-13))];
julia> Z=[transpose([LinRange(m-r+1e-2,m+r-1e-2,250)' m-r.+2*r*(exp.(im*LinRange(0,pi,250)')/2 .+.5)]); Z[1:500]];
julia> z,fz,ω,err=svAAA(nep.nep2,Z);
julia> err[end]
4.589829024476524e-14
```
# References
- S. Güttel, G. M. Negri Porzio, and F. Tisseur. Robust approximations of
  nonlinear eigenvalue problems. 2020
- P. Lietaert, J. Pérez, B. Vandereycken, and K. Meerbergen. Automatic rational
  approximation and linearization of nonlinear eigenvalue problems. IMA journal
  of numerical analysis, 2018.
- Y. Nakatsukasa, O. Sète, and L. N. Trefethen. The AAA algorithm for rational
  approximation. SIAM journal on scientific computing, 40(3):A1494-A1522, 2018.
"""
svAAA(nep::AbstractSPMF, Z::AbstractArray{T}; params...) where T<:Number =
    svAAA(ComplexF64, nep, Z; params...)
function svAAA(
        ::Type{T},
        nep::AbstractSPMF,
        Z::AbstractArray{T};
        mmax::Int = 100,
        tol::Real = eps(real(T))*1e3,
        cleanup::Bool = true,
        tol_cln::Real = min(eps(real(T)),tol),
        return_details::Bool = false,
        logger = 0,
        weighted::Bool = false) where T<:Number

    @parse_logger_param!(logger)

    fv = get_fv(nep)
    Z = Z[:]
    Z = deleteat!(Z,isinf.(Z) .|| isnan.(Z))
    M = length(Z)
    s = length(fv)
    F = [fi(z) for z in Z, fi in fv]

    # AAA-variant-dependent parameters
    if weighted
        Av = get_Av(nep)
        n = size(nep,1)
        rng = copy(Random.default_rng())
        Random.seed!(1)
        u = randn(T,n)
        copy!(Random.default_rng(),rng)
        u = u/norm(u)
        uj = Matrix{T}(undef,n,s)
        for j = 1:s
            mul!((@view uj[:,j]),Av[j],u)
        end
        beta = 0
        for i = 1:M
            beta = max(beta, norm(uj*( @view F[i,:] )))
        end

        scaleF = transpose(map(i -> norm(Av[i]), 1:s))
        F = F.*scaleF
        scaleF = 1 ./scaleF
        if cleanup  # For normalization of residues
            maxF = maximum(abs.(F),dims=1)
        end
    else
        scaleF = maximum(abs.(F),dims=1)
        F = F./scaleF
    end

    # Initializations of arrays
    err = Vector{real(T)}(undef,0)
    z = Vector{T}(undef,0)
    ω = Vector{T}(undef,0)
    fz = Matrix{T}(undef,mmax,s)
    ind = Vector{Int}(undef,0)

    H = UpperTriangular(Matrix{T}(undef,mmax,mmax))
    S = UpperTriangular(Matrix{T}(undef,mmax,mmax))

    Q = Matrix{T}(undef,M*s,min(50,mmax))
    C = Matrix{T}(undef,M,min(50,mmax))

    vmat = Matrix{T}(undef,M,s)

    res = Matrix{real(T)}(undef,M,s)
    N = Matrix{T}(undef,M,s)
    D = Vector{T}(undef,M)
    R = Matrix{T}(undef,M,s)
    R .= mean(F,dims=1) # Initial value for residuals

    for m = 1:mmax
        # Finding new support point
        res .= abs.(F.-R)
        maxres, idx = findmax(res,dims=1)
        loc = idx[argmax(maxres)]
        locz = loc[1]
        push!(err, ( weighted ? sum(maxres)/beta : res[loc] ))
        push_info!(logger, 2, "svAAA iteration $(m-1): Error = $(err[m])")
        if err[m] <= tol # Convergence check
            fz = scaleF.*( @view fz[1:m-1,:] )
            break
        end
        push!(z,Z[locz])
        push!(ind,locz)
        fz[m,:] = @view F[locz,:]

        # Find new column of Loewner matrix
        #   expansion of Cauchy matrix
        if size(C,2) < m
            resized = zeros(T,M,min(2*(m-1),mmax))
            resized[1:end, 1:(m-1)] = C
            C = resized
        end
        C[:,m] .= one(T)./(Z.-z[m])
        C[ind,m] .= 0
        vmat .= @views ( C[:,m].*(F .- transpose(fz[m,:])) )
        v = vec(vmat)

        # (Partial) update of QR-decomposition
        q = @views (Q[locz:M:end,1:m-1]*S[1:m-1,1:m-1] )
        ee = I - q'*q
        chol = cholesky!(ee)
        Si = chol.U
        H[1:m-1,1:m-1] = Si*( @view H[1:m-1,1:m-1] )
        S[1:m-1,1:m-1] = ( @view S[1:m-1,1:m-1] )/Si
        S[m,1:m-1] .= 0
        S[1:m-1,m] .= 0
        S[m,m] = 1
        Q[locz:M:end,1:m-1] .= 0

        # Gram-Schmidt orthogonalization
        nv = norm(v)
        mul!(( @view H[1:m-1,m] ), ( @view Q[:,1:m-1] )', v)
        H[1:m-1,m] = @views ( S[1:m-1,1:m-1]'*H[1:m-1,m] )
        HH = @views ( S[1:m-1,1:m-1]*H[1:m-1,m] )
        mul!(v, ( @view Q[:,1:m-1] ), HH, -one(T), one(T))
        H[m,m] = norm(v)
        #   DGKS-reorthogonalization
        ii = 0
        while (ii < 3) && (real(H[m,m]) < 1/sqrt(2)*nv)
            mul!(HH, ( @view Q[:,1:m-1])', v)
            HH = @views ( S[1:m-1,1:m-1]'*HH )
            H[1:m-1,m] .= ( @view H[1:m-1,m] ) .+ HH
            HH = @views ( S[1:m-1,1:m-1]*HH )
            mul!(v, ( @view Q[:,1:m-1] ), HH, -one(T), one(T))
            nv = real(H[m,m])
            H[m,m] = norm(v)
            ii = ii+1
        end
        #   expansion of Q-matrix
        if size(Q,2) < m
            resized = zeros(T, M*s, min(2*(m-1),mmax))
            resized[1:end, 1:(m-1)] = Q
            Q = resized
        end
        Q[:,m] .= v./H[m,m]

        # Solve least-squares using SVD
        SVD = svd!(H[1:m,1:m])
        ω = SVD.V[:,end]

        # Update values in sample set
        mul!(N, ( @view C[:,1:m] ), ( ω.*( @view fz[1:m,:] ) ))
        mul!(D, ( @view C[:,1:m] ), ω)
        R .= N./D
        R[ind,:] .= @view F[ind,:]

        # Detection/removal of spurious poles
        if cleanup && m>1
            # Compute poles via generalized eigenvalue problem
            B = Matrix{T}(I,m+1,m+1)
            B[1,1] = 0 ;
            E = [ 0 transpose(ω) ; ones(m,1) diagm(z)]
            eig = eigen!(E,B)
            pol = eig.values
            pol = deleteat!(pol,isnan.(pol))

            # Compute residues via discretized Cauchy integral
            dz = 1e-5*[im, -1, -im, 1]
            pp = pol .+ transpose(dz)
            rvals = reval(reshape(pp,(m-1)*4), z, (@view fz[1:m,:]), ω)
            rsd = Matrix{T}(undef,length(pol),s)
            for i = 1:s
                rsd[:,i] = reshape((@view rvals[:,i]), (length(pol), 4) )*dz./4
            end

            # Check criterium for spurious poles
            if weighted
                maxRsd = maximum(abs, rsd./maxF, dims=2)
            else # If SV-AAA is used, functions are already of the same scale
                maxRsd = maximum(abs, rsd, dims=2)
            end
            idx = findall(maxRsd.< tol_cln)
            nb_sp = length(idx)

            # Remove support points closest to spurious poles (if any)
            if nb_sp > 0
                loc_sp = Vector{Int}(undef,nb_sp)
                ind_sp = Vector{Int}(undef,nb_sp)
                for j = 1:nb_sp
                    azp = abs.(z .- pol[idx[j]])
                    locj = argmin(azp)
                    loc_sp[j] = locj
                    ind_sp[j] = ind[locj]
                    deleteat!(z,locj)
                    deleteat!(ind,locj)
                end
                loc_z = setdiff(1:m,loc_sp)
                ind_Z = setdiff(1:M,ind)

                # Find new weights using SVD of Loewner matrix
                C[ind_sp,loc_z] .= @views ( 1 ./ (Z[ind_sp] .- transpose(z)) ) # These values were zero
                Cview = @view C[ind_Z, loc_z]
                L = Matrix{T}(undef,(M-m+nb_sp)*s,m-nb_sp)
                for j=1:s
                    L[1+(j-1)*(M-m+nb_sp):j*(M-m+nb_sp),:] .= Cview.*(F[ind_Z,j] .- transpose(fz[loc_z,j]))
                end
                SVD = svd!(L)
                ω = SVD.V[:,end]

                # Update values in sample set
                N .= @views ( C[:,loc_z]*(ω.*fz[loc_z,:]) )
                D .= ( @view C[:,loc_z] )*ω
                R .= N./D
                R[ind,:] .= @view F[ind,:]

                # Compute last error and exit for-loop
                res = abs.(F.-R)
                maxres, idx = findmax(res,dims=1)
                loc = idx[argmax(maxres)]
                push!(err, ( weighted ? sum(maxres)/beta : res[loc] ))
                if nb_sp == 1
                    push_info!(logger, "svAAA: 1 Froissart doublet detected (and removed). Final error = $(err[m+1])")
                else
                    push_info!(logger, "svAAA: $nb_sp Froissart doublets detected (and removed). Final error = $(err[m+1])")
                end
                fz = scaleF.*( @view fz[loc_z, :] )
                break
            end
        end

        # End of iterations
        if m == mmax
            res .= abs.(F.-R)
            maxres, idx = findmax(res,dims=1)
            loc = idx[argmax(maxres)]
            push!(err, ( weighted ? sum(maxres)/beta : res[loc] ))
            if err[m+1] <= tol
                push_info!(logger, 2, "svAAA iteration $mmax: Error = $(err[mmax+1])")
            else
                push_info!(logger, "svAAA: Rational approximation not converged after $mmax iterations. Final error = $(err[mmax+1])")
            end
            fz = scaleF.*( @view fz[1:m,:] )
        end
    end

    # Remove support points with zero weight
    idx_z = findall(ω.==0)
    if length(idx_z) > 0
        deleteat!(z,idx_z)
        fz = fz[setdiff(1:end,idx_z),:]
        deleteat!(ω,idx_z)
    end

    # Poles, residues and zeros
    if return_details
        pol, rsd, zer = get_prz(z, fz, ω)
    else
        pol = Vector{T}() ; rsd = Vector{T}() ; zer = Vector{T}()
    end

    return z, fz, ω, err, pol, rsd, zer

end

# Evaluate rational interpolant(s) (defined by z, fz and ω) in λ
function reval(λ, z::AbstractVector{T}, fz::AbstractArray{T}, ω::AbstractVector{T}) where T<:Number
    C = one(T)./(λ.-transpose(z))
    r = ( C*(ω.*fz) )./(C*ω)

    # Deal with input inf: r(inf) = lim r(zz) = sum(ω.*f) / sum(ω):
    r[isinf.(λ),:] = repeat( sum(ω.*fz,dims=1)./sum(ω) , inner=( sum(isinf.(λ)) , 1 ) )

    # Deal with NaN:
    # i1, i2 = findnz(isnan.(r)) # findnz is not included for BitArray (Plots adds own version for dense arrays) -> findall(iszero,...)
    keysnz = findall(!iszero,isnan.(r))
    i1 = [k[1] for k in keysnz]
    i2 = [k[2] for k in keysnz]

    for j = 1:size(i1,1)
        if ( isnan(λ[i1[j]]) || ~any((λ[i1[j]]) .== z) )
            # r(NaN) = NaN is fine.
            # The second case may happen if r(zv(ii)) = 0/0 at some point.
        else
            # Clean up values NaN = inf/inf at support points.
            r[i1[j],i2[j]] = fz[findfirst(isone,λ[i1[j]] .== z),i2[j]];
        end
    end
    return r
end

# Return poles, residues and zeros of rational interpolant(s) defined by (z, fz and ω)
function get_prz(z::AbstractVector{T}, fz::AbstractArray{T}, ω::AbstractVector{T}) where T<:Number
    m, s = size(fz)

    # Compute poles via generalized eigenvalue problem
    B = Matrix{T}(I,m+1,m+1)
    B[1,1] = 0 ;
    E = [ 0 transpose(ω) ; ones(m,1) diagm(z)]
    eig = eigen(E,B)
    pol = eig.values
    pol = pol[(~).(isnan.(pol))]

    # Compute residues via discretized Cauchy integral
    dz = 1e-5*[im, -1, -im, 1]
    pp = pol .+ transpose(dz)
    rvals = reval(reshape(pp,length(pol)*4), z, (@view fz[1:m,:]), ω)
    rsd = Matrix{T}(undef,length(pol),s)
    for i = 1:s
        rsd[:,i] = reshape(rvals[:,i], (length(pol), 4) )*dz./4
    end

    # Compute zeros via generalized eigenvalue problem
    zer = Matrix{T}(undef,m+1,s)
    for i = 1:s
        E[1,2:end] .= ω.*( @view fz[:,i] )
        eig = eigen(E,B)
        zer[:,i] = eig.values
    end

    return pol, rsd, zer
end

# Struct compiling optional output of AAAeigs
struct AAASolutionDetails{T<:Number}
    m_appr::Int
    zfω::AbstractArray
    prz::AbstractArray
    err_appr::AbstractVector
    Lam::AbstractArray{T}
    Res::AbstractArray
    conv_it::Int
end

# Empty constructor
AAASolutionDetails{T}() where T<:Number = AAASolutionDetails(
    0,
    Vector{AbstractMatrix{T}}(),
    Vector{AbstractMatrix{T}}(),
    Vector{Real}(),
    Matrix{T}(undef,0,0),
    Matrix{Real}(undef,0,0),
    0)
