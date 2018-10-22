using Pkg, Dates, Printf;
Pkg.activate("."); # Assumes . is the NonlinearEigenproblems directory
Pkg.instantiate();

hpc_script="hpc1.jl"; # Script to be run
include(hpc_script);

# Create a directory for the profiling data:
datetime=mkdir(Dates.format(now(),"yyyy-mm-dd-HH:MM:SS"));
profile_data_dir="prf_"*datetime*"_";
profile_data_dir=profile_data_dir*replace(hpc_script,".jl" => "")
mkdir(profile_data_dir)

println("profile data dir: ",profile_data_dir);

# Copy the HPC-script to the profiling directory
s = read(joinpath(@__DIR__(),hpc_script), String)
open(joinpath(profile_data_dir,"hpc_script.jl"), "w") do f
    write(f, s)
end


# Profile the code
using Profile
Profile.init(n=10^7,delay=0.005)
# From documentation: n is the total number of instruction pointers you can store, with a default value of 10^6. If your typical backtrace is 20 instruction pointers, then you can collect 50000 backtraces, which suggests a statistical uncertainty of less than 1%. This may be good enough for most applications. The default setting is delay = 0.001. Of...

# n=10^7 delay=0.005 => 2 hours runtime

t1 = time_ns()
run_hpc(preprun=true) # preparation run without profiling first to compile, etc
t2 = time_ns()
@profile run_hpc()
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
