module OpenWaterEvap

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export OWEVAP_DATA,
       OWEVAP_HARGREAVES_1985,
       OWEVAP_HARGREAVES,
       OWEVAP_PENMAN_SIMPLE,
       OWEVAP_PENMAN_MONTEITH,
       OPEN_WATER_EVAP,
       OPEN_WATER_RIPARIAN,
       OPEN_WATER_UWFS

"""Open-water evaporation provided directly by forcing."""
function OWEVAP_DATA(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    evaporation_forcing::Number=first(@variables evaporation_forcing),
    waterstorage::Number=first(@variables waterstorage),
    flux_name::Symbol=:owevap_data,
)
    @hydroflux flux_name openwater_evaporation ~ clamp(max(0.0, evaporation_forcing), 0.0, max(0.0, waterstorage))
end

"""Classic Hargreaves (1985) open-water evaporation."""
function OWEVAP_HARGREAVES_1985(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    mean_temp::Number=first(@variables mean_temp),
    max_temp::Number=first(@variables max_temp),
    min_temp::Number=first(@variables min_temp),
    extraterrestrial_radiation::Number=first(@variables extraterrestrial_radiation),
    waterstorage::Number=first(@variables waterstorage),
    flux_name::Symbol=:owevap_hargreaves_1985,
)
    evap = 0.0023 * max(0.0, extraterrestrial_radiation) * (mean_temp + 17.8) * sqrt(max(max_temp - min_temp, 0.0))
    @hydroflux flux_name openwater_evaporation ~ clamp(evap, 0.0, max(0.0, waterstorage))
end

"""Generalized Hargreaves open-water evaporation."""
function OWEVAP_HARGREAVES(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    mean_temp::Number=first(@variables mean_temp),
    max_temp::Number=first(@variables max_temp),
    min_temp::Number=first(@variables min_temp),
    extraterrestrial_radiation::Number=first(@variables extraterrestrial_radiation),
    waterstorage::Number=first(@variables waterstorage),
    kh::Number=first(@parameters kh [description = "Hargreaves coefficient", bounds = (0, 1), unit = "-"]),
    ch::Number=first(@parameters ch [description = "Hargreaves temperature offset", bounds = (-50, 50), unit = "degC"]),
    flux_name::Symbol=:owevap_hargreaves,
)
    evap = kh * max(0.0, extraterrestrial_radiation) * (mean_temp + ch) * sqrt(max(max_temp - min_temp, 0.0))
    @hydroflux flux_name openwater_evaporation ~ clamp(evap, 0.0, max(0.0, waterstorage))
end

"""Simple Penman open-water evaporation."""
function OWEVAP_PENMAN_SIMPLE(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    net_radiation::Number=first(@variables net_radiation),
    slope_sat_vapour_curve::Number=first(@variables slope_sat_vapour_curve),
    psychrometric_constant::Number=first(@variables psychrometric_constant),
    latent_heat::Number=first(@variables latent_heat),
    wind_function::Number=first(@variables wind_function),
    vapour_pressure_deficit::Number=first(@variables vapour_pressure_deficit),
    waterstorage::Number=first(@variables waterstorage),
    flux_name::Symbol=:owevap_penman_simple,
)
    denom = max(slope_sat_vapour_curve + psychrometric_constant, 1.0e-12)
    evap = (slope_sat_vapour_curve / denom) * (net_radiation / max(latent_heat, 1.0e-12)) +
           (psychrometric_constant / denom) * wind_function * vapour_pressure_deficit
    @hydroflux flux_name openwater_evaporation ~ clamp(evap, 0.0, max(0.0, waterstorage))
end

"""Penman-Monteith open-water evaporation."""
function OWEVAP_PENMAN_MONTEITH(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    net_radiation::Number=first(@variables net_radiation),
    ground_heat_flux::Number=first(@variables ground_heat_flux),
    slope_sat_vapour_curve::Number=first(@variables slope_sat_vapour_curve),
    psychrometric_constant::Number=first(@variables psychrometric_constant),
    latent_heat::Number=first(@variables latent_heat),
    air_density::Number=first(@variables air_density),
    specific_heat_air::Number=first(@variables specific_heat_air),
    vapour_pressure_deficit::Number=first(@variables vapour_pressure_deficit),
    aerodynamic_resistance::Number=first(@variables aerodynamic_resistance),
    waterstorage::Number=first(@variables waterstorage),
    flux_name::Symbol=:owevap_penman_monteith,
)
    numerator = slope_sat_vapour_curve * (net_radiation - ground_heat_flux) +
                air_density * specific_heat_air * vapour_pressure_deficit / max(aerodynamic_resistance, 1.0e-12)
    denominator = max(latent_heat, 1.0e-12) * max(slope_sat_vapour_curve + psychrometric_constant, 1.0e-12)
    @hydroflux flux_name openwater_evaporation ~ clamp(numerator / denominator, 0.0, max(0.0, waterstorage))
end

"""Basic open-water evaporation as a scaled PET over open water."""
function OPEN_WATER_EVAP(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    pet_open_water::Number=first(@variables pet_open_water),
    waterstorage::Number=first(@variables waterstorage),
    scale_coeff::Number=first(@parameters scale_coeff [description = "Open-water PET scaling coefficient", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:open_water_evap,
)
    @hydroflux flux_name openwater_evaporation ~ clamp(scale_coeff * max(0.0, pet_open_water), 0.0, max(0.0, waterstorage))
end

"""Riparian open-water evaporation scaled by a saturated-area fraction."""
function OPEN_WATER_RIPARIAN(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    pet_open_water::Number=first(@variables pet_open_water),
    sat_fraction::Number=first(@variables sat_fraction),
    waterstorage::Number=first(@variables waterstorage),
    scale_coeff::Number=first(@parameters scale_coeff [description = "Open-water PET scaling coefficient", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:open_water_riparian,
)
    @hydroflux flux_name openwater_evaporation ~ clamp(scale_coeff * max(0.0, sat_fraction) * max(0.0, pet_open_water), 0.0, max(0.0, waterstorage))
end

"""UWFS-style wetland/depression open-water evaporation scaled by depression fraction."""
function OPEN_WATER_UWFS(;
    openwater_evaporation::Number=first(@variables openwater_evaporation),
    pet_open_water::Number=first(@variables pet_open_water),
    depression_fraction::Number=first(@variables depression_fraction),
    waterstorage::Number=first(@variables waterstorage),
    scale_coeff::Number=first(@parameters scale_coeff [description = "Open-water PET scaling coefficient", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:open_water_uwfs,
)
    @hydroflux flux_name openwater_evaporation ~ clamp(scale_coeff * max(0.0, depression_fraction) * max(0.0, pet_open_water), 0.0, max(0.0, waterstorage))
end

end
