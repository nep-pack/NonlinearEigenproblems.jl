#!/bin/bash
# Submit benchmark created by run_benchmark.jl to GitHub. Intended to be run by Travis.
set -ev
branch=`git rev-parse --abbrev-ref HEAD`; if [ $branch = 'HEAD' ]; then branch=`git describe --all --exact-match 2>/dev/null | sed 's=.*/=='`; fi
commit=`git rev-parse --short HEAD`
benchmark_file_name=benchmark-$branch-$TRAVIS_BUILD_NUMBER-$commit.json
git clone --depth=10 --branch=master https://github.com/nep-pack/julia-ci.git julia-ci-benchmark
mv benchmark.json julia-ci-benchmark/$benchmark_file_name
cd julia-ci-benchmark
git config user.email benchmark
git config user.name benchmark
git add $benchmark_file_name
git commit -m "benchmark $commit"
git push "https://$GIT_TOKEN@github.com/nep-pack/julia-ci.git" master > /dev/null 2>&1
