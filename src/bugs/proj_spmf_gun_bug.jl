workspace()
push!(LOAD_PATH, string(@__DIR__, "/.."))
push!(LOAD_PATH, string(@__DIR__, "/../gallery_extra"))

using NEPCore
using NEPTypes
using Gallery
using GalleryNLEVP

nep_org = nep_gallery(NLEVP_NEP,"gun")
nep = nlevp_make_native(nep_org)

proj_nep = create_proj_NEP(nep)
V = rand(9956, 10);

set_projectmatrices!(proj_nep, V, V)
