#!/bin/bash
julia --color=yes -e 'using Pkg; Pkg.activate("test"); include("test/benchmark/print_benchmarks.jl")' -- "$@"
