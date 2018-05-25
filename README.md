
This package is a first alpha version of methods for nonlinear eigenvalue problems in julia.

Core developers (alphabetical): Elias Jarlebring, Giampaolo Mele, Emil Ringh, Parikshit Upadhyaya

Links below are on KTH we will migrate them soon: 

Documentation page: https://gitr.sys.kth.se/pages/nep-pack/nep-pack-alpha/


* Checklist for first public version: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/issues/26

* NEP-page style "guide": https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/Style-guidelines-and-notes

* GIT-usage: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki

* NEP-methods: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/NEP-methods

* NEP-applications: https://gitr.sys.kth.se/nep-pack/nep-pack-alpha/wiki/Applications



julia> if is_windows(); Pkg.add("WinRPM"); using WinRPM; WinRPM.install("gcc"); end
julia> include(joinpath(dirname(JULIA_HOME), "share", "julia", "build_sysimg.jl"))
julia> build_sysimg(force=true)
