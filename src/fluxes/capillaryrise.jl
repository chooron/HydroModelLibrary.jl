module CapillaryRise

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export CAPRISE_HBV,
       CRISE_HBV

"""HBV-style capillary rise from a lower store into the soil store."""
function CAPRISE_HBV(;
    capillary_rise::Number=first(@variables capillary_rise),
    lower_storage::Number=first(@variables lower_storage),
    soil_storage::Number=first(@variables soil_storage),
    max_soil_storage::Number=first(@parameters max_soil_storage [description = "Maximum soil storage", bounds = (0, 5000), unit = "mm"]),
    capillary_coeff::Number=first(@parameters capillary_coeff [description = "Maximum capillary rise coefficient", bounds = (0, 100), unit = "mm/d"]),
    flux_name::Symbol=:caprise_hbv,
)
    demand_factor = max(1 - max(0.0, soil_storage) / max(max_soil_storage, 1.0e-12), 0.0)
    @hydroflux flux_name capillary_rise ~ min(max(0.0, lower_storage), capillary_coeff * demand_factor)
end

"""Alias matching the Raven Chapter 3 process name."""
function CRISE_HBV(; kwargs...)
    params = (; kwargs...)
    return CAPRISE_HBV(; params..., flux_name=get(params, :flux_name, :crise_hbv))
end

end
