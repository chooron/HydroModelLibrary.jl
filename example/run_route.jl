include("../src/HydroModelLibrary.jl")
include("../src/routes/gamma_uh.jl")
using HydroModels
using ComponentArrays
@variables flow flow_routed
@parameters α β

uh = GammaUnitHydro.GammaHydrograph([flow], [flow_routed], α, β)

input = reshape(Float64[1, 3, 5, 7, 8, 22, 34, 45, 32, 21, 11, 4], 1, :)
params = ComponentVector(params=(α=0.5, β=0.5))
output = uh(input, params)