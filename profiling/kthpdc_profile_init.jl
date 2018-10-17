using Pkg, Dates;
Pkg.activate("."); # Assumes . is the NonlinearEigenproblems directory
Pkg.instantiate();


hpc_script="hpc2.jl"; # Script to be run

# Create a directory for the profiling data:
datetime=mkdir(Dates.format(now(),"yyyy-mm-dd-HH:MM:SS"));
profile_data_dir="prf_"*datetime*"_";
profile_data_dir=profile_data_dir*replace(hpc_script,".jl" => "")
mkdir(profile_data_dir)

# Profile the code
using Profile
include(hpc_script) # always run without profiling first to compile, etc
@profile include(hpc_script);

# Extract profiling data to variables
li, lidict = Profile.retrieve()

# Store it to disk (non-compatible version)
using Serialization

open(joinpath(profile_data_dir,"profile_li.data"), "w") do f; serialize(f, li); end
open(joinpath(profile_data_dir,"profile_lidict.data"), "w") do f; serialize(f, lidict); end
