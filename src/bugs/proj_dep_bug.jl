workspace()
push!(LOAD_PATH, string(@__DIR__, "/.."))
using NEPCore
using NEPTypes
using Gallery

nep = nep_gallery("dep0",100)
proj_nep = create_proj_NEP(nep)

V = rand(100,5);

set_projectmatrices!(proj_nep, V, V)

compute_Mder(proj_nep, 1+0.5im , 1)
