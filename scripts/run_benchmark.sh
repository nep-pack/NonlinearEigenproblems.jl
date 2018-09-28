#!/bin/bash
julia -O3 --project=test --color=yes test/benchmark/run_benchmark.jl -- "$@"
