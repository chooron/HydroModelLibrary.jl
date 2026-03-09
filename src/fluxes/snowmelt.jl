module SnowMelt

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export SNOWMELT_DEGREE_DAY, SNOWMELT_FROM_POTENTIAL

"""Actual snowmelt limited by snow storage and supplied potential melt."""
function SNOWMELT_FROM_POTENTIAL(;
    snowmelt::Number=first(@variables snowmelt),
    snowstorage::Number=first(@variables snowstorage),
    potential_melt::Number=first(@variables potential_melt),
    flux_name::Symbol=:snowmelt_from_potential,
)
    @hydroflux flux_name snowmelt ~ min(max(0.0, snowstorage), max(0.0, potential_melt))
end

"""Actual snowmelt from a degree-day formulation capped by snow storage."""
function SNOWMELT_DEGREE_DAY(;
    snowmelt::Number=first(@variables snowmelt),
    snowstorage::Number=first(@variables snowstorage),
    temp::Number=first(@variables temp),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (-10, 10), unit = "degC"]),
    ma::Number=first(@parameters ma [description = "Melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    flux_name::Symbol=:snowmelt_degree_day,
)
    @hydroflux flux_name snowmelt ~ min(max(0.0, snowstorage), ma * max(0.0, temp - tf))
end

end
