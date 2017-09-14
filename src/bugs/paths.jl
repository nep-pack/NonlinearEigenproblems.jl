# TEST: Test different ways to produce paths.
# When files/modules load eachother within the package it is important to be precise about paths.
# Not a bug in Julia, more of an illustration of things that can produce bugs in our package.

workspace()

println("")
println("=================================================")
println("||                 TESTING PATHS!              ||")
println("=================================================")
println("")
println("The first line will change depending on where your current directory is (where you started Julia).")
println("   pwd() = ", pwd())
println("The second line is where the current file is located (a macro).")
println("   @__DIR__ = ", @__DIR__)

println("")
println("The first changes, like after calling cd(\"..\") as we do now")

bname = basename(pwd())
cd("..")

println("   pwd() = ", pwd())
println("But the second is unaffected of this call")
println("   @__DIR__ = ", @__DIR__)

cd(bname) #Go back to where we started
