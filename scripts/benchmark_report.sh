#!/bin/bash
julia --project=test --color=yes test/benchmark/benchmark_report.jl -- "$@"
