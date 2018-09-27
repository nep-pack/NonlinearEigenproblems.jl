#!/bin/bash
julia --color=yes -e 'using Pkg; Pkg.activate("test"); include("test/runtests.jl")' -- "$@"
