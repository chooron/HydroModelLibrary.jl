module CanopyEvap

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export CANEVAP_ALL,
       CANEVAP_LINEAR,
       CANEVAP_PRMS,
       CANEVP_MAXIMUM,
       CANEVP_RUTTER

"""Canopy evaporation limited only by canopy storage and PET."""
function CANEVAP_ALL(;
    canopy_evaporation::Number=first(@variables canopy_evaporation),
    canopy_storage::Number=first(@variables canopy_storage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:canevap_all,
)
    @hydroflux flux_name canopy_evaporation ~ clamp(max(0.0, potential_evaporation * pet_corr), 0.0, max(0.0, canopy_storage))
end

"""Linear canopy evaporation reduction with relative canopy storage."""
function CANEVAP_LINEAR(;
    canopy_evaporation::Number=first(@variables canopy_evaporation),
    canopy_storage::Number=first(@variables canopy_storage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    max_canopy_storage::Number=first(@parameters max_canopy_storage [description = "Maximum canopy storage", bounds = (0, 5000), unit = "mm"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:canevap_linear,
)
    storage_ratio = clamp(max(0.0, canopy_storage) / max(max_canopy_storage, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name canopy_evaporation ~ clamp(max(0.0, potential_evaporation * pet_corr) * storage_ratio, 0.0, max(0.0, canopy_storage))
end

"""PRMS-style canopy evaporation using a vegetation-area fraction."""
function CANEVAP_PRMS(;
    canopy_evaporation::Number=first(@variables canopy_evaporation),
    canopy_storage::Number=first(@variables canopy_storage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    vegetation_fraction::Number=first(@parameters vegetation_fraction [description = "Vegetated area fraction", bounds = (0, 1), unit = "-"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:canevap_prms,
)
    @hydroflux flux_name canopy_evaporation ~ clamp(max(0.0, vegetation_fraction * potential_evaporation * pet_corr), 0.0, max(0.0, canopy_storage))
end

"""Maximum canopy evaporation from PET, canopy fraction, and snow-free fraction."""
function CANEVP_MAXIMUM(;
    canopy_evaporation::Number=first(@variables canopy_evaporation),
    potential_evaporation::Number=first(@variables potential_evaporation),
    canopy_fraction::Number=first(@variables canopy_fraction),
    snow_fraction::Number=first(@variables snow_fraction),
    canopy_storage::Number=first(@variables canopy_storage),
    flux_name::Symbol=:canevp_maximum,
)
    @hydroflux flux_name canopy_evaporation ~ clamp(max(0.0, potential_evaporation) * clamp(canopy_fraction, 0.0, 1.0) * (1 - clamp(snow_fraction, 0.0, 1.0)), 0.0, max(0.0, canopy_storage))
end

"""Rutter canopy evaporation scaled by canopy storage relative to canopy capacity."""
function CANEVP_RUTTER(;
    canopy_evaporation::Number=first(@variables canopy_evaporation),
    potential_evaporation::Number=first(@variables potential_evaporation),
    canopy_fraction::Number=first(@variables canopy_fraction),
    trunk_fraction::Number=first(@variables trunk_fraction),
    canopy_storage::Number=first(@variables canopy_storage),
    canopy_capacity::Number=first(@parameters canopy_capacity [description = "Canopy storage capacity", bounds = (1.0e-6, 5000), unit = "mm"]),
    flux_name::Symbol=:canevp_rutter,
)
    storage_ratio = clamp(max(0.0, canopy_storage) / max(canopy_capacity, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name canopy_evaporation ~ clamp(max(0.0, potential_evaporation) * clamp(canopy_fraction, 0.0, 1.0) * (1 - clamp(trunk_fraction, 0.0, 1.0)) * storage_ratio, 0.0, max(0.0, canopy_storage))
end

end
