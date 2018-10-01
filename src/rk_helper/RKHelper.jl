"""
Rational Krylov helper module. Contains functionality needed for methods
like NLEIGS and CORK.
"""
module RKHelper

export inpolygon, discretizepolygon,
       lejabagby, ratnewtoncoeffs, ratnewtoncoeffsm,
       LinSolverCache, solve

include("inpolygon.jl")
include("discretizepolygon.jl")
include("lejabagby.jl")
include("ratnewtoncoeffs.jl")
include("linsolvercache.jl")

end
