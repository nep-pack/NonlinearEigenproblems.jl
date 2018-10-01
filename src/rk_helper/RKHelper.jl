"""
Rational Krylov helper module. Contains functionality needed for methods
like NLEIGS and CORK.
"""
module RKHelper

export inpolygon, discretizepolygon,
       lejabagby, scgendivdiffs, ratnewtoncoeffs, ratnewtoncoeffsm,
       LinSolverCache, solve

include("inpolygon.jl")
include("discretizepolygon.jl")
include("rk_utils.jl")
include("linsolvercache.jl")

end
