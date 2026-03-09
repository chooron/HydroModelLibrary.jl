module CanopyDrip

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export CANDRIP_EXCESS,
       CANDRIP_PRMS,
       CANDRIP_LINEAR,
       CANDRIP_SLOWDRAIN

"""Canopy drip / throughfall as storage above canopy capacity."""
function CANDRIP_EXCESS(;
    canopy_drip::Number=first(@variables canopy_drip),
    canopy_storage::Number=first(@variables canopy_storage),
    max_canopy_storage::Number=first(@parameters max_canopy_storage [description = "Maximum canopy storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:candrip_excess,
)
    @hydroflux flux_name canopy_drip ~ max(max(0.0, canopy_storage) - max_canopy_storage, 0.0)
end

"""PRMS-style throughfall driven by incoming intercepted precipitation once canopy storage exceeds capacity."""
function CANDRIP_PRMS(;
    canopy_drip::Number=first(@variables canopy_drip),
    canopy_storage::Number=first(@variables canopy_storage),
    intercepted_precipitation::Number=first(@variables intercepted_precipitation),
    max_canopy_storage::Number=first(@parameters max_canopy_storage [description = "Maximum canopy storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:candrip_prms,
)
    drip_trigger = ifelse(canopy_storage > max_canopy_storage, 1.0, 0.0)
    @hydroflux flux_name canopy_drip ~ min(max(0.0, canopy_storage), drip_trigger * max(0.0, intercepted_precipitation))
end

"""Linear canopy drainage above a storage threshold."""
function CANDRIP_LINEAR(;
    canopy_drip::Number=first(@variables canopy_drip),
    canopy_storage::Number=first(@variables canopy_storage),
    drip_coeff::Number=first(@parameters drip_coeff [description = "Canopy drainage coefficient", bounds = (0, 10), unit = "d-1"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Canopy drainage threshold", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:candrip_linear,
)
    @hydroflux flux_name canopy_drip ~ clamp(drip_coeff * max(max(0.0, canopy_storage) - storage_threshold, 0.0), 0.0, max(0.0, canopy_storage))
end

"""Slow canopy drainage proportional to relative canopy storage."""
function CANDRIP_SLOWDRAIN(;
    canopy_drip::Number=first(@variables canopy_drip),
    canopy_storage::Number=first(@variables canopy_storage),
    canopy_capacity::Number=first(@parameters canopy_capacity [description = "Canopy capacity", bounds = (1.0e-6, 5000), unit = "mm"]),
    drip_coeff::Number=first(@parameters drip_coeff [description = "Slow-drain coefficient", bounds = (0, 100), unit = "mm/d"]),
    flux_name::Symbol=:candrip_slowdrain,
)
    @hydroflux flux_name canopy_drip ~ clamp(drip_coeff * max(0.0, canopy_storage) / max(canopy_capacity, 1.0e-12), 0.0, max(0.0, canopy_storage))
end

end
