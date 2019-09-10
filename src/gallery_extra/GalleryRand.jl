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

export MSWS_RNG
export gen_rng_int
export gen_rng_float
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
    reset_rng!(seed::UInt128) = reset_rng!(GLOBAL_MSWS_RNG, seed);
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
    gen_rng_float(lower::Float64, upper::Float64) = gen_rng_float(GLOBAL_MSWS_RNG, lower::Float64, upper::Float64)
    function gen_rng_float(rng::MSWS_RNG, lower::Float64, upper::Float64)
        return upper - (upper-lower)*gen_rng_float(rng);
    end

end
