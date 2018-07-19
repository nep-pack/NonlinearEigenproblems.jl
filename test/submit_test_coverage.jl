###########################################
# Submits test coverage for relevant code #
###########################################
using Coverage

src_folder = normpath(string(@__DIR__) * "/../src")

# Add folders and jl files below if you want to exclude them from the test
# coverage reports. Folders are excluded recursively. All names are case
# insensitive, and relative to the src/ folder.
excluded_folders_and_files = map(f -> uppercase(src_folder * "/" * f), [
    "bugs/",
    "extra_tests/",
    "tmp/",
    "trash/",
    ])

coverage = process_folder(src_folder)

filter!(f -> !any(x -> startswith(uppercase(f.filename), x), excluded_folders_and_files), coverage)

Coveralls.submit(coverage)
Codecov.submit(coverage)
