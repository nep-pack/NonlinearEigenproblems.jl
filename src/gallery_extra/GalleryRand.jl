"""
    GalleryRand

A module which handles generation of pseudo random matrices.
Allows stability over releases.
B. Widynski, Middle Square Weyl Sequence RNG, arXiv 1704.00358

# Reference implementation:
    mutable struct MSWS_RNG_32
        x::UInt64
        w::UInt64
        s::UInt64
        function MSWS_RNG_32()
            #this = new(0,0,0xb5ad4eceda1ce2a9)
            this = new(0,0,0x0000000100000001)
            return this
        end
    end

    function gen_rng_int(rng::MSWS_RNG_32)
        rng.x *= rng.x
        rng.x += (rng.w += rng.s)
        rng.x = (rng.x>>32) | (rng.x<<32)
        return UInt32(rng.x & 0xFFFFFFFF)
    end
"""
module GalleryRand

using LinearAlgebra
using SparseArrays

export MSWS_RNG
export gen_rng_int
export gen_rng_float
export gen_rng_mat
export gen_rng_spmat


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

    function gen_rng_mat(rng::MSWS_RNG, n::Int64, m::Int64)
        A = zeros(Float64,n,m)
        for c = 1:m
            for r = 1:n
                A[r,c] = 1-2*gen_rng_float(rng)
            end
        end
        return A
    end

    function gen_rng_spmat(rng::MSWS_RNG, n::Int64, m::Int64, p::Real)
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

end
