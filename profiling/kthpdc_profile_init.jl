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
t1 = time_ns()
include(hpc_script) # always run without profiling first to compile, etc
t2 = time_ns()
@profile include(hpc_script);
t3 = time_ns()

# Extract profiling data to variables
li, lidict = Profile.retrieve()

# Store it to disk (non-compatible version)
using Serialization

open(joinpath(profile_data_dir,"profile_li.data"), "w") do f; serialize(f, li); end
open(joinpath(profile_data_dir,"profile_lidict.data"), "w") do f; serialize(f, lidict); end

include("metadata_serialization.jl")
file_name = joinpath(profile_data_dir, "metadata.json")
total_time = (t3 - t1) / 1e9
profile_time = (t3 - t2) / 1e9
# TODO: pass any other relevant custom fields to the method below, such as the
# hpc_script, or a manually entered comment
save_metadata(file_name, total_time, profile_time)
