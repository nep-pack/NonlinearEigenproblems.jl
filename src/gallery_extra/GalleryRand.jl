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
            #this = new(0,0,0xb5ad4eceda1ce2a9);
            this = new(0,0,0x0000000100000001);
            return this
        end
    end

    function gen_rng_int(rng::MSWS_RNG_32)
        rng.x *= rng.x;
        rng.x += (rng.w += rng.s);
        rng.x = (rng.x>>32) | (rng.x<<32);
        return UInt32(rng.x & 0xFFFFFFFF)
    end
"""
module GalleryRand


using LinearAlgebra
using SparseArrays


export MSWS_RNG
export gen_rng
export gen_rng_int
export gen_rng_float
export gen_rng_complex
export gen_rng_mat
export gen_rng_spmat
export reset_rng!


    mutable struct MSWS_RNG
        x::UInt128
        w::UInt128
        s::UInt128
        function MSWS_RNG(seed::UInt128 = zero(UInt128))
            this = new();
            reset_rng!(this, seed);
            return this
        end
    end


    reset_rng!() = reset_rng!(GLOBAL_MSWS_RNG, zero(UInt128));
    reset_rng!(rng::MSWS_RNG, seed::Integer) = reset_rng!(rng, unsigned(Int128(seed)));
    reset_rng!(seed::Integer) = reset_rng!(GLOBAL_MSWS_RNG, seed);
    function reset_rng!(rng::MSWS_RNG, seed::UInt128)
        base = 0x9ef09a97ac0f9ecaef01c4f2db0958c9; # Found by a random draw
        rng.s = (seed << 1) + base; # Left-shift to make sure seed i even, then add to base which is odd (https://github.com/tidwall/weyl/issues/1)
        rng.x = 0x1de568e1a1ca1b593cbf13f7407cf43e;
        rng.w = 0xd4ac5c288559e14a5fafc1b7df9f9e0e;
        return
    end



    const GLOBAL_MSWS_RNG = MSWS_RNG();



    function gen_rng_int(rng::MSWS_RNG = GLOBAL_MSWS_RNG)
        rng.x *= rng.x;
        rng.x += (rng.w += rng.s);
        rng.x = (rng.x>>64) | (rng.x<<64);
        return UInt64(rng.x & typemax(UInt64))
    end


    gen_rng_int(lower::Integer, upper::Integer) = gen_rng_int(GLOBAL_MSWS_RNG, lower, upper)
    function gen_rng_int(rng::MSWS_RNG, lower::Integer, upper::Integer)
        if lower >= upper
            error("Lower limit cannot be larger than or equal to upper limit.")
        elseif lower < 0
            return signed(gen_rng_int(rng, UInt64(upper-lower))) + lower
        else
            return gen_rng_int(rng, UInt64(upper-lower)) + lower
        end
    end


    gen_rng_int(upper::Integer) = gen_rng_int(GLOBAL_MSWS_RNG, UInt64(unsigned(upper)));
    gen_rng_int(rng::MSWS_RNG, upper::Integer) = gen_rng_int(rng, UInt64(unsigned(upper)));
    gen_rng_int(upper::Unsigned) = gen_rng_int(GLOBAL_MSWS_RNG, upper);
    gen_rng_int(rng::MSWS_RNG, upper::Unsigned) = gen_rng_int(rng, UInt64(upper));
    gen_rng_int(upper::UInt64) = gen_rng_int(GLOBAL_MSWS_RNG, lower, upper)
    function gen_rng_int(rng::MSWS_RNG, upper::UInt64)
        return mod(gen_rng_int(rng), upper+1) #Include the upper limit, but wrong if upper=typemax(UInt64)
    end


    function gen_rng_float(rng::MSWS_RNG = GLOBAL_MSWS_RNG)
        return Float64(gen_rng_int(rng)/typemax(UInt64));
    end


    gen_rng_float(lower::Real, upper::Real) = gen_rng_float(GLOBAL_MSWS_RNG, Float64(lower), Float64(upper))
    gen_rng_float(lower::Float64, upper::Float64) = gen_rng_float(GLOBAL_MSWS_RNG, lower, upper)
    gen_rng_float(rng::MSWS_RNG, lower::Real, upper::Real) = gen_rng_float(rng, Float64(lower), Float64(upper))
    function gen_rng_float(rng::MSWS_RNG, lower::Float64, upper::Float64)
        return upper - (upper-lower)*gen_rng_float(rng);
    end


    gen_rng_complex(lower::Real, upper::Complex) = gen_rng_complex(GLOBAL_MSWS_RNG, ComplexF64(lower), ComplexF64(upper))
    gen_rng_complex(lower::Complex, upper::Real) = gen_rng_complex(GLOBAL_MSWS_RNG, ComplexF64(lower), ComplexF64(upper))
    gen_rng_complex(lower::Complex, upper::Complex) = gen_rng_complex(GLOBAL_MSWS_RNG, ComplexF64(lower), ComplexF64(upper))
    gen_rng_complex(lower::ComplexF64, upper::ComplexF64) = gen_rng_complex(GLOBAL_MSWS_RNG, lower, upper)
    gen_rng_complex(rng::MSWS_RNG, lower::Real, upper::Complex) = gen_rng_complex(rng, ComplexF64(lower), ComplexF64(upper))
    gen_rng_complex(rng::MSWS_RNG, lower::Complex, upper::Real) = gen_rng_complex(rng, ComplexF64(lower), ComplexF64(upper))
    gen_rng_complex(rng::MSWS_RNG, lower::Complex, upper::Complex) = gen_rng_complex(rng, ComplexF64(lower), ComplexF64(upper))
    function gen_rng_complex(rng::MSWS_RNG, lower::ComplexF64, upper::ComplexF64)
        # Defines a rectange in complex plane from which we draw real and imag uniformly
        return gen_rng_float(rng, real(lower), real(upper)) + one(ComplexF64)im*gen_rng_float(rng, imag(lower), imag(upper))
    end


    gen_rng(lower::Tl,upper::Tu) where {Tl,Tu} = gen_rng(Tl,GLOBAL_MSWS_RNG,lower,upper)
    gen_rng(::Type{T},lower::Tl,upper::Tu) where {T,Tl,Tu} = gen_rng(T,GLOBAL_MSWS_RNG,lower,upper)
    gen_rng(::Type{Int64},rng::MSWS_RNG,lower::Int64,upper::Int64) = gen_rng_int(rng,lower,upper)
    gen_rng(::Type{Float64},rng::MSWS_RNG,lower::Float64,upper::Float64) = gen_rng_float(rng,lower,upper)
    gen_rng(::Type{ComplexF64},rng::MSWS_RNG,lower::ComplexF64,upper::ComplexF64) = gen_rng_complex(rng,lower,upper)
    function gen_rng(::Type{T},rng::MSWS_RNG,lower::T,upper::T) where T
        error("This should not happen. Not implemented for type ", T)
    end
    function gen_rng(::Type{T},rng::MSWS_RNG,lower::Tl,upper::Tu) where {T,Tl,Tu}
        TT = promote_type(Tl,Tu)
        if (TT<:Complex); TT = ComplexF64; end;
        gen_rng(TT, rng, TT(lower), TT(upper))
    end


    gen_rng_mat(n::Int64, m::Int64, lower, upper) = gen_rng_mat(GLOBAL_MSWS_RNG, n, m, lower, upper)
    gen_rng_mat(rng::MSWS_RNG, n::Int64, lower, upper) = gen_rng_mat(rng, n, n, lower, upper)
    gen_rng_mat(n::Int64, lower, upper) = gen_rng_mat(GLOBAL_MSWS_RNG, n, n, lower, upper)
    function gen_rng_mat(rng::MSWS_RNG, n::Int64, m::Int64, lower::Tl, upper::Tu) where Tl where Tu
        T = promote_type(Tl,Tu)
        if (T<:Complex); T = ComplexF64; end;
        gen_rng_mat(rng, n, m, T(lower), T(upper))
    end
    function gen_rng_mat(rng::MSWS_RNG, n::Int64, m::Int64, lower::T, upper::T) where T <: Union{Int64,Float64,ComplexF64}
        A = zeros(T,n,m)
        for c = 1:m
            for r = 1:n
                A[r,c] = gen_rng(T,rng,lower,upper)
            end
        end
        return A
    end


    gen_rng_spmat(n::Int64, m::Int64, p::Real, lower, upper) = gen_rng_spmat(GLOBAL_MSWS_RNG, n, m, p, lower, upper)
    gen_rng_spmat(rng::MSWS_RNG, n::Int64, p::Real, lower, upper) = gen_rng_spmat(rng, n, n, p, lower, upper)
    gen_rng_spmat(n::Int64, p::Real, lower, upper) = gen_rng_spmat(GLOBAL_MSWS_RNG, n, n, p, lower, upper)
    function gen_rng_spmat(rng::MSWS_RNG, n::Int64, m::Int64, p::Real, lower::Tl, upper::Tu) where Tl where Tu
        T = promote_type(Tl,Tu)
        if (T<:Complex); T = ComplexF64; end;
        gen_rng_spmat(rng, n, m, p, T(lower), T(upper))
    end
    function gen_rng_spmat(rng::MSWS_RNG, n::Int64, m::Int64, p::Real, lower::T, upper::T) where T <: Union{Int64,Float64,ComplexF64}
        nonzeros = round(p*m*n)
        dict = Dict{Tuple{Int64,Int64},T}()
        for i = 1:nonzeros
            r = gen_rng_int(rng,n-1)+1
            c = gen_rng_int(rng,m-1)+1
            dict[r,c] = gen_rng(T,rng,lower,upper)
        end
        idxes = collect(keys(dict))
        nonzeros = length(idxes)
        r = zeros(Int64,nonzeros)
        c = zeros(Int64,nonzeros)
        vals = zeros(T,nonzeros)
        for i = 1:nonzeros
            r[i] = idxes[i][1]
            c[i] = idxes[i][2]
            vals[i] = dict[idxes[i]]
        end
        return sparse(r,c,vals,n,m)
    end


end
