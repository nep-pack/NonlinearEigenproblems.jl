#!/bin/bash
julia --project=test --color=yes test/runtests.jl -- "$@"
