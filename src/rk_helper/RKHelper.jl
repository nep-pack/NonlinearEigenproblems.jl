"""
Rational Krylov helper module. Contains functionality needed for methods
like NLEIGS and CORK.
"""
module RKHelper

include("rk_nep.jl")
include("rk_utils.jl")
include("linearize.jl")

include("inpolygon.jl")
include("discretizepolygon.jl")

include("linsolvercache.jl")

end
