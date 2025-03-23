module HydroModelLibrary

using HydroModels
using ComponentArrays
using ModelingToolkit
using Lux
using Reexport
using NamedTupleTools

fluxes_names = [
    "baseflow", "capillary", "evaporation", "exchange", "infiltration", "interception", "interflow",
    "melt", "normalize", "percolation", "rainfall", "recharge", "refreeze", "saturation", "snowfall"
]
bucket_names = ["cemaneige", "exphydro", "gr4j", "hymod", "simhyd", "m50", "m100", "dplHBV", "hbv"]

map(fluxes_names) do name
    include("fluxes/$(name).jl")
end

map(bucket_names) do name
    include("buckets/$(name).jl")
end





include("marrmot/collie1.jl")
export Collie1

export Cemaneige, ExpHydro, GR4J, HyMOD, SIMHYD, M50, M100, HBV
end
