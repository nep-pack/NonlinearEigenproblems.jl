#!/bin/bash
julia --color=yes -e 'using Pkg; Pkg.activate("."); include("test/benchmark/benchmark_report.jl")' "$@"
