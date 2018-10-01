#!/bin/bash
julia --project=test --color=yes test/benchmark/print_benchmarks.jl -- "$@"
