module PrecipInterception

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export PRECIP_ICEPT_LAI,
       PRECIP_ICEPT_EXPLAI

"""LAI-based linear canopy interception partition for rainfall and snowfall."""
function PRECIP_ICEPT_LAI(;
    intercepted_rain::Number=first(@variables intercepted_rain),
    intercepted_snow::Number=first(@variables intercepted_snow),
    rainfall::Number=first(@variables rainfall),
    snowfall::Number=first(@variables snowfall),
    lai::Number=first(@variables lai),
    sai::Number=first(@variables sai),
    alpha_rain::Number=first(@parameters alpha_rain [description = "Rain interception coefficient", bounds = (0, 10), unit = "-"]),
    alpha_snow::Number=first(@parameters alpha_snow [description = "Snow interception coefficient", bounds = (0, 10), unit = "-"]),
    rain_flux_name::Symbol=:precip_icept_lai_rain,
    snow_flux_name::Symbol=:precip_icept_lai_snow,
)
    cover_index = max(0.0, lai + sai)
    rain_fraction = clamp(alpha_rain * cover_index, 0.0, 1.0)
    snow_fraction = clamp(alpha_snow * cover_index, 0.0, 1.0)
    rain_flux = @hydroflux rain_flux_name intercepted_rain ~ max(0.0, rainfall) * rain_fraction
    snow_flux = @hydroflux snow_flux_name intercepted_snow ~ max(0.0, snowfall) * snow_fraction
    return (rain_flux, snow_flux)
end

"""Exponential-LAI canopy interception partition for rainfall and snowfall."""
function PRECIP_ICEPT_EXPLAI(;
    intercepted_rain::Number=first(@variables intercepted_rain),
    intercepted_snow::Number=first(@variables intercepted_snow),
    rainfall::Number=first(@variables rainfall),
    snowfall::Number=first(@variables snowfall),
    lai::Number=first(@variables lai),
    sai::Number=first(@variables sai),
    rain_flux_name::Symbol=:precip_icept_explai_rain,
    snow_flux_name::Symbol=:precip_icept_explai_snow,
)
    cover_index = max(0.0, lai + sai)
    interception_fraction = clamp(1 - exp(-0.5 * cover_index), 0.0, 1.0)
    rain_flux = @hydroflux rain_flux_name intercepted_rain ~ max(0.0, rainfall) * interception_fraction
    snow_flux = @hydroflux snow_flux_name intercepted_snow ~ max(0.0, snowfall) * interception_fraction
    return (rain_flux, snow_flux)
end

end
