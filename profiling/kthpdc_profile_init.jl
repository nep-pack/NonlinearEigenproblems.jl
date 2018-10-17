using Pkg;
Pkg.activate("."); # Assumes . is the NonlinearEigenproblems directory
Pkg.instantiate();


# Profile some code
using Profile
hpc_script="hpc3.jl"; # Script to be run
include(hpc_script) # always run without profiling first to compile, etc
@profile include(hpc_script)

# Extract profiling data to variables
li, lidict = Profile.retrieve()

# Store it to disk (non-compatible version)
using Serialization
open("profile_li.data", "w") do f; serialize(f, li); end
open("profile_lidict.data", "w") do f; serialize(f, lidict); end
