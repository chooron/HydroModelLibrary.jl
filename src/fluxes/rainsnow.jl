module RainSnow

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export RAINSNOW_DATA, RAINSNOW_DINGMAN, RAINSNOW_HBV, RAINSNOW_UBCWM

"""Rain/snow partition provided directly by forcing inputs."""
function RAINSNOW_DATA(;
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    snowfall_forcing::Number=first(@variables snowfall_forcing),
    rainfall_forcing::Number=first(@variables rainfall_forcing),
    snow_flux_name::Symbol=:rainsnow_data_snow,
    rain_flux_name::Symbol=:rainsnow_data_rain,
)
    snow_flux = @hydroflux snow_flux_name snowfall ~ max(0.0, snowfall_forcing)
    rain_flux = @hydroflux rain_flux_name rainfall ~ max(0.0, rainfall_forcing)
    return (snow_flux, rain_flux)
end

"""Dingman rain/snow partition."""
function RAINSNOW_DINGMAN(;
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    precipitation::Number=first(@variables precipitation),
    temp::Number=first(@variables temp),
    threshold_temp::Number=first(@parameters threshold_temp [description = "Dingman transition temperature", bounds = (-10, 10), unit = "degC"]),
    snow_flux_name::Symbol=:rainsnow_dingman_snow,
    rain_flux_name::Symbol=:rainsnow_dingman_rain,
)
    cold_term = 1.0 - 0.5 * exp(-2.2 * max(threshold_temp - temp, 0.0)^1.3)
    warm_term = 0.5 * exp(-2.2 * max(temp - threshold_temp, 0.0)^1.3)
    snow_fraction = clamp(ifelse(temp <= threshold_temp, cold_term, warm_term), 0.0, 1.0)
    snow_flux = @hydroflux snow_flux_name snowfall ~ max(0.0, precipitation) * snow_fraction
    rain_flux = @hydroflux rain_flux_name rainfall ~ max(0.0, precipitation) * (1 - snow_fraction)
    return (snow_flux, rain_flux)
end

"""HBV rain/snow partition."""
function RAINSNOW_HBV(;
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    precipitation::Number=first(@variables precipitation),
    temp::Number=first(@variables temp),
    tt::Number=first(@parameters tt [description = "Snowfall temperature threshold", bounds = (-10, 10), unit = "degC"]),
    tti::Number=first(@parameters tti [description = "Rain/snow transition interval", bounds = (1.0e-6, 20), unit = "degC"]),
    snow_flux_name::Symbol=:rainsnow_hbv_snow,
    rain_flux_name::Symbol=:rainsnow_hbv_rain,
)
    snow_fraction = clamp((tt + tti / 2 - temp) / max(tti, 1.0e-12), 0.0, 1.0)
    snow_flux = @hydroflux snow_flux_name snowfall ~ max(0.0, precipitation) * snow_fraction
    rain_flux = @hydroflux rain_flux_name rainfall ~ max(0.0, precipitation) * (1 - snow_fraction)
    return (snow_flux, rain_flux)
end

"""UBCWM-style partition using mean temperature and daily temperature range."""
function RAINSNOW_UBCWM(;
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    precipitation::Number=first(@variables precipitation),
    temp::Number=first(@variables temp),
    delta_temp::Number=first(@variables delta_temp),
    snow_temp::Number=first(@parameters snow_temp [description = "All-snow threshold temperature", bounds = (-10, 10), unit = "degC"]),
    rain_temp::Number=first(@parameters rain_temp [description = "All-rain threshold temperature", bounds = (-10, 10), unit = "degC"]),
    snow_flux_name::Symbol=:rainsnow_ubcwm_snow,
    rain_flux_name::Symbol=:rainsnow_ubcwm_rain,
)
    tmax = temp + 0.5 * max(0.0, delta_temp)
    tmin = temp - 0.5 * max(0.0, delta_temp)
    snow_fraction = ifelse(
        tmax <= snow_temp,
        1.0,
        ifelse(
            tmin >= rain_temp,
            0.0,
            clamp((rain_temp - temp) / max(rain_temp - snow_temp, 1.0e-12), 0.0, 1.0),
        ),
    )
    snow_flux = @hydroflux snow_flux_name snowfall ~ max(0.0, precipitation) * snow_fraction
    rain_flux = @hydroflux rain_flux_name rainfall ~ max(0.0, precipitation) * (1 - snow_fraction)
    return (snow_flux, rain_flux)
end

end