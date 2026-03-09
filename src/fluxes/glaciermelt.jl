module GlacierMelt

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export GLACIERMELT_DEGREE_DAY,
       GLACIERMELT_FROM_POTENTIAL,
       GLACIERMELT_GSMSOCONT

_snow_free_factor(snow_storage::Number, snow_shield_scale::Number) =
    1 / (1 + exp((max(0.0, snow_storage) - 0.05 * max(snow_shield_scale, 1.0e-12)) / (0.01 * max(snow_shield_scale, 1.0e-12))))

"""Actual glacier melt from a standard degree-day formulation capped by available glacier storage."""
function GLACIERMELT_DEGREE_DAY(;
    glacier_melt::Number=first(@variables glacier_melt),
    glacier_storage::Number=first(@variables glacier_storage),
    temp::Number=first(@variables temp),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (-10, 10), unit = "degC"]),
    ma::Number=first(@parameters ma [description = "Glacier melt factor", bounds = (0, 50), unit = "mm/(d*degC)"]),
    flux_name::Symbol=:glaciermelt_degree_day,
)
    melt = ma * max(0.0, temp - tf)
    @hydroflux flux_name glacier_melt ~ min(max(0.0, glacier_storage), max(0.0, melt))
end

"""Actual glacier melt limited by supplied potential melt and available glacier storage."""
function GLACIERMELT_FROM_POTENTIAL(;
    glacier_melt::Number=first(@variables glacier_melt),
    glacier_storage::Number=first(@variables glacier_storage),
    potential_melt::Number=first(@variables potential_melt),
    flux_name::Symbol=:glaciermelt_from_potential,
)
    @hydroflux flux_name glacier_melt ~ min(max(0.0, glacier_storage), max(0.0, potential_melt))
end

"""GSMSOCONT-style glacier melt suppressed by overlying snow cover."""
function GLACIERMELT_GSMSOCONT(;
    glacier_melt::Number=first(@variables glacier_melt),
    glacier_storage::Number=first(@variables glacier_storage),
    snow_storage::Number=first(@variables snow_storage),
    temp::Number=first(@variables temp),
    tm::Number=first(@parameters tm [description = "Glacier melt temperature threshold", bounds = (-10, 10), unit = "degC"]),
    aice::Number=first(@parameters aice [description = "Degree-day factor for ice melt", bounds = (0, 50), unit = "mm/(d*degC)"]),
    snow_shield_scale::Number=first(@parameters snow_shield_scale [description = "Snow storage scale suppressing ice melt", bounds = (1.0e-6, 5000), unit = "mm"]),
    flux_name::Symbol=:glaciermelt_gsmsocont,
)
    melt = _snow_free_factor(snow_storage, snow_shield_scale) * aice * max(0.0, temp - tm)
    @hydroflux flux_name glacier_melt ~ min(max(0.0, glacier_storage), max(0.0, melt))
end

end
