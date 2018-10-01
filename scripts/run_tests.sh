#!/bin/bash
julia --color=yes -e 'using Pkg; Pkg.activate("."); include("test/runtests.jl")' -- "$@"
