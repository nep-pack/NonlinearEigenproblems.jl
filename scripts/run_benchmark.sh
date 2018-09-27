#!/bin/bash
julia -O3 --color=yes -e 'using Pkg; Pkg.activate("test"); include("test/benchmark/run_benchmark.jl")' -- "$@"
