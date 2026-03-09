module SnowAlbedo

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export SNOALB_UBC,
       SNOALB_CRHM_ESSERY,
       SNOALB_BAKER

function SNOALB_UBC(;
    snow_albedo_change::Number=first(@variables snow_albedo_change),
    snow_albedo::Number=first(@variables snow_albedo),
    snow_input::Number=first(@variables snow_input),
    albedo_max::Number=first(@parameters albedo_max [description = "Maximum snow albedo", bounds = (0, 1), unit = "-"]),
    decay_alpha::Number=first(@parameters decay_alpha [description = "UBC albedo decay parameter", bounds = (0, 1), unit = "-"]),
    retain_k::Number=first(@parameters retain_k [description = "UBC albedo retention parameter", bounds = (0, 1), unit = "-"]),
    snow_albedo_scale::Number=first(@parameters snow_albedo_scale [description = "Snowfall albedo renewal scale", bounds = (1.0e-6, 5000), unit = "mm"]),
    dt::Number=first(@variables dt [description = "Time step in days", bounds = (1.0e-6, 365)]),
    flux_name::Symbol=:snoalb_ubc,
)
    decay_term = -decay_alpha * (1 - retain_k) / max(dt, 1.0e-12)
    renewal_term = (albedo_max - snow_albedo) / max(dt, 1.0e-12) * min(max(0.0, snow_input) / max(snow_albedo_scale, 1.0e-12), 1.0)
    @hydroflux flux_name snow_albedo_change ~ decay_term + renewal_term
end

function SNOALB_CRHM_ESSERY(;
    snow_albedo_change::Number=first(@variables snow_albedo_change),
    snow_albedo::Number=first(@variables snow_albedo),
    snow_temp::Number=first(@variables snow_temp),
    snowfall::Number=first(@variables snowfall),
    albedo_min::Number=first(@parameters albedo_min [description = "Minimum snow albedo", bounds = (0, 1), unit = "-"]),
    albedo_max::Number=first(@parameters albedo_max [description = "Maximum snow albedo", bounds = (0, 1), unit = "-"]),
    albedo_break::Number=first(@parameters albedo_break [description = "Albedo breakpoint", bounds = (0, 1), unit = "-"]),
    beta_cold::Number=first(@parameters beta_cold [description = "Cold-snow decay parameter", bounds = (0, 1), unit = "-"]),
    beta_warm::Number=first(@parameters beta_warm [description = "Warm-snow decay parameter", bounds = (0, 1), unit = "-"]),
    snowfall_threshold::Number=first(@parameters snowfall_threshold [description = "Snowfall threshold for albedo reset", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:snoalb_crhm_essery,
)
    decay = ifelse(
        snow_temp < 0,
        -beta_cold,
        ifelse(snow_albedo < albedo_break, -beta_warm * (snow_albedo - albedo_min), 0.0),
    )
    renewal = (albedo_max - snow_albedo) * min(max(0.0, snowfall) / max(snowfall_threshold, 1.0e-12), 1.0)
    @hydroflux flux_name snow_albedo_change ~ decay + renewal
end

function SNOALB_BAKER(;
    snow_albedo::Number=first(@variables snow_albedo),
    snow_age::Number=first(@variables snow_age),
    flux_name::Symbol=:snoalb_baker,
)
    @hydroflux flux_name snow_albedo ~ 0.9 - 0.0473 * max(snow_age, 0.0)^0.1
end

end
