module SnowFreeze

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export FREEZE_DEGREE_DAY, FREEZE_HMETS

"""Degree-day refreezing."""
function FREEZE_DEGREE_DAY(;
    snowfreeze::Number=first(@variables snowfreeze),
    liquidwater::Number=first(@variables liquidwater),
    temp::Number=first(@variables temp),
    kf::Number=first(@parameters kf [description = "Degree-day factor for freezing", bounds = (0, 10), unit = "mm/(d*degC)"]),
    tf::Number=first(@parameters tf [description = "Freezing temperature", bounds = (-10, 10), unit = "degC"]),
    flux_name::Symbol=:freeze_degree_day,
)
    @hydroflux flux_name snowfreeze ~ min(max(0.0, liquidwater), kf * max(0.0, tf - temp))
end

"""HMETS refreezing using an overnight temperature proxy."""
function FREEZE_HMETS(;
    snowfreeze::Number=first(@variables snowfreeze),
    liquidwater::Number=first(@variables liquidwater),
    mean_temp::Number=first(@variables mean_temp),
    min_temp::Number=first(@variables min_temp),
    kf::Number=first(@parameters kf [description = "HMETS degree-day factor for refreezing", bounds = (0, 10), unit = "mm/(d*degC^f)"]),
    tbf::Number=first(@parameters tbf [description = "Refreezing temperature threshold", bounds = (-10, 10), unit = "degC"]),
    freezing_n::Number=first(@parameters freezing_n [description = "HMETS refreezing exponent", bounds = (0, 5), unit = "-"]),
    flux_name::Symbol=:freeze_hmets,
)
    overnight_temp = (mean_temp + min_temp) / 2
    @hydroflux flux_name snowfreeze ~ min(max(0.0, liquidwater), kf * max(0.0, tbf - overnight_temp)^freezing_n)
end

end
