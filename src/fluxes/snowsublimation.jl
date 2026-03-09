module SnowSublimation

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export SUBLIM_KUZMIN,
       SUBLIM_CENTRAL_SIERRA,
       SUBLIM_BULK_AERO

function SUBLIM_KUZMIN(;
    snow_sublimation::Number=first(@variables snow_sublimation),
    wind_speed::Number=first(@variables wind_speed),
    sat_vapour_pressure::Number=first(@variables sat_vapour_pressure),
    vapour_pressure::Number=first(@variables vapour_pressure),
    flux_name::Symbol=:sublim_kuzmin,
)
    @hydroflux flux_name snow_sublimation ~ max(0.0, (0.18 + 0.098 * max(0.0, wind_speed)) * max(sat_vapour_pressure - vapour_pressure, 0.0))
end

function SUBLIM_CENTRAL_SIERRA(;
    snow_sublimation::Number=first(@variables snow_sublimation),
    wind_speed::Number=first(@variables wind_speed),
    sat_vapour_pressure::Number=first(@variables sat_vapour_pressure),
    vapour_pressure::Number=first(@variables vapour_pressure),
    hw::Number=first(@parameters hw [description = "Surface roughness height", bounds = (1.0e-6, 1.0e6), unit = "m"]),
    hv::Number=first(@parameters hv [description = "Vegetation roughness height", bounds = (1.0e-6, 1.0e6), unit = "m"]),
    flux_name::Symbol=:sublim_central_sierra,
)
    coeff = 0.0063 * (max(hw * hv, 1.0e-12))^(-1 / 6)
    @hydroflux flux_name snow_sublimation ~ max(0.0, coeff * max(sat_vapour_pressure - vapour_pressure, 0.0) * max(0.0, wind_speed))
end

function SUBLIM_BULK_AERO(;
    snow_sublimation::Number=first(@variables snow_sublimation),
    air_density::Number=first(@variables air_density),
    latent_heat_sublimation::Number=first(@variables latent_heat_sublimation),
    wind_speed::Number=first(@variables wind_speed),
    sat_vapour_pressure_snow::Number=first(@variables sat_vapour_pressure_snow),
    air_vapour_pressure::Number=first(@variables air_vapour_pressure),
    air_pressure::Number=first(@variables air_pressure),
    water_density::Number=first(@variables water_density),
    zref::Number=first(@parameters zref [description = "Reference height", bounds = (1.0e-6, 1.0e6), unit = "m"]),
    z0::Number=first(@parameters z0 [description = "Momentum roughness length", bounds = (1.0e-6, 1.0e6), unit = "m"]),
    z0e::Number=first(@parameters z0e [description = "Vapour roughness length", bounds = (1.0e-6, 1.0e6), unit = "m"]),
    flux_name::Symbol=:sublim_bulk_aero,
)
    kappa = 0.41
    ce = kappa^2 / max(log(zref / max(z0, 1.0e-12)) * log(zref / max(z0e, 1.0e-12)), 1.0e-12)
    qe = air_density * latent_heat_sublimation * ce * max(0.0, wind_speed) * 0.622 * max(sat_vapour_pressure_snow - air_vapour_pressure, 0.0) / max(air_pressure, 1.0e-12)
    @hydroflux flux_name snow_sublimation ~ max(0.0, qe / max(latent_heat_sublimation * water_density, 1.0e-12))
end

end
