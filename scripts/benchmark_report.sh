#!/bin/bash
julia --color=yes -e 'using Pkg; Pkg.activate("test"); include("test/benchmark/benchmark_report.jl")' -- "$@"
