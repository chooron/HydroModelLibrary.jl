module HydroModelLibrary

using HydroModels

# fluxes_names = [
#     "baseflow", "capillary", "evaporation", "exchange", "infiltration", "interception", "interflow",
#     "melt", "normalize", "percolation", "rainfall", "recharge", "refreeze", "saturation", "snowfall"
# ]
# bucket_names = ["cemaneige", "exphydro", "gr4j", "hymod", "simhyd", "m50", "m100", "dplHBV", "hbv"]

# map(fluxes_names) do name
#     include("fluxes/$(name).jl")
# end

# map(bucket_names) do name
#     include("buckets/$(name).jl")
# end
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

include("models/alpine1.jl")
include("models/collie2.jl")
include("models/gr4j.jl")
include("models/susannah2.jl")
include("models/wetland.jl")
export alpine1
end
