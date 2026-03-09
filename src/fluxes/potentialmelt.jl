module PotentialMelt

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export POTMELT_DATA,
       POTMELT_DEGREE_DAY,
       POTMELT_DD_FREEZE,
       POTMELT_HBV,
       POTMELT_HBV_ROS,
       POTMELT_HMETS,
       POTMELT_RESTRICTED

"""Potential snowmelt prescribed directly by forcing."""
function POTMELT_DATA(;
    potential_melt::Number=first(@variables potential_melt),
    melt_forcing::Number=first(@variables melt_forcing),
    flux_name::Symbol=:potmelt_data,
)
    @hydroflux flux_name potential_melt ~ max(0.0, melt_forcing)
end

"""Potential snowmelt using a standard degree-day formulation."""
function POTMELT_DEGREE_DAY(;
    potential_melt::Number=first(@variables potential_melt),
    temp::Number=first(@variables temp),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (-10, 10), unit = "degC"]),
    ma::Number=first(@parameters ma [description = "Melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    flux_name::Symbol=:potmelt_degree_day,
)
    @hydroflux flux_name potential_melt ~ ma * max(0.0, temp - tf)
end

"""Signed degree-day melt/freezing potential."""
function POTMELT_DD_FREEZE(;
    potential_melt::Number=first(@variables potential_melt),
    temp::Number=first(@variables temp),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (-10, 10), unit = "degC"]),
    ma::Number=first(@parameters ma [description = "Degree-day factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    flux_name::Symbol=:potmelt_dd_freeze,
)
    @hydroflux flux_name potential_melt ~ ma * (temp - tf)
end

"""HBV potential melt with an optional multiplicative correction factor."""
function POTMELT_HBV(;
    potential_melt::Number=first(@variables potential_melt),
    temp::Number=first(@variables temp),
    melt_corr::Number=first(@variables melt_corr [description = "Optional HBV melt correction factor", bounds = (0, 10)]),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (-10, 10), unit = "degC"]),
    ma::Number=first(@parameters ma [description = "Melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    flux_name::Symbol=:potmelt_hbv,
)
    @hydroflux flux_name potential_melt ~ ma * max(0.0, temp - tf) * max(0.0, melt_corr)
end

"""HBV rain-on-snow potential melt with an additive rainfall correction."""
function POTMELT_HBV_ROS(;
    potential_melt::Number=first(@variables potential_melt),
    temp::Number=first(@variables temp),
    rainfall::Number=first(@variables rainfall),
    melt_corr::Number=first(@variables melt_corr [description = "Optional HBV melt correction factor", bounds = (0, 10)]),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (-10, 10), unit = "degC"]),
    ma::Number=first(@parameters ma [description = "Melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    ros_coeff::Number=first(@parameters ros_coeff [description = "Rain-on-snow melt coefficient", bounds = (0, 100), unit = "mm-1"]),
    flux_name::Symbol=:potmelt_hbv_ros,
)
    thermal_term = max(0.0, temp - tf)
    @hydroflux flux_name potential_melt ~ max(0.0, ma * max(0.0, melt_corr) * thermal_term + ros_coeff * max(0.0, rainfall) * thermal_term)
end

"""HMETS potential melt with a cumulative-melt-dependent degree-day factor."""
function POTMELT_HMETS(;
    potential_melt::Number=first(@variables potential_melt),
    temp::Number=first(@variables temp),
    cumulmelt::Number=first(@variables cumulmelt),
    mamax::Number=first(@parameters mamax [description = "Maximum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    mamin::Number=first(@parameters mamin [description = "Minimum melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    malpha::Number=first(@parameters malpha [description = "Cumulative melt factor coefficient", bounds = (0, 1), unit = "-"]),
    tbm::Number=first(@parameters tbm [description = "Melt temperature threshold", bounds = (-10, 10), unit = "degC"]),
    flux_name::Symbol=:potmelt_hmets,
)
    melt_factor = min(mamax, mamin * (1 + malpha * max(0.0, cumulmelt)))
    @hydroflux flux_name potential_melt ~ melt_factor * max(0.0, temp - tbm)
end

"""Degree-day melt restricted by currently available snow and snowfall."""
function POTMELT_RESTRICTED(;
    potential_melt::Number=first(@variables potential_melt),
    temp::Number=first(@variables temp),
    snowstorage::Number=first(@variables snowstorage),
    snowfall::Number=first(@variables snowfall),
    melt_corr::Number=first(@variables melt_corr [description = "Potential melt correction factor", bounds = (0, 10)]),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (-10, 10), unit = "degC"]),
    ma::Number=first(@parameters ma [description = "Melt factor", bounds = (0, 20), unit = "mm/(d*degC)"]),
    flux_name::Symbol=:potmelt_restricted,
)
    degree_day_melt = max(0.0, melt_corr) * ma * max(0.0, temp - tf)
    available_snow = max(0.0, snowstorage + snowfall)
    @hydroflux flux_name potential_melt ~ min(available_snow, degree_day_melt)
end

end
