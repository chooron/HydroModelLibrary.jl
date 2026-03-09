module GlacierRelease

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export GRELEASE_LINEAR_STORAGE,
       GRELEASE_LINEAR_ANALYTIC,
       GRELEASE_HBV_EC

function GRELEASE_LINEAR_STORAGE(;
    glacier_release::Number=first(@variables glacier_release),
    glacier_storage::Number=first(@variables glacier_storage),
    release_coeff::Number=first(@parameters release_coeff [description = "Linear glacier-release coefficient", bounds = (0, 10), unit = "d-1"]),
    flux_name::Symbol=:grelease_linear_storage,
)
    @hydroflux flux_name glacier_release ~ max(0.0, release_coeff * glacier_storage)
end

function GRELEASE_LINEAR_ANALYTIC(;
    glacier_release::Number=first(@variables glacier_release),
    glacier_storage::Number=first(@variables glacier_storage),
    release_coeff::Number=first(@parameters release_coeff [description = "Linear glacier-release coefficient", bounds = (0, 10), unit = "d-1"]),
    dt::Number=first(@variables dt [description = "Time step in days", bounds = (1.0e-6, 365)]),
    flux_name::Symbol=:grelease_linear_analytic,
)
    @hydroflux flux_name glacier_release ~ max(0.0, glacier_storage) / max(dt, 1.0e-12) * (1 - exp(-release_coeff * dt))
end

function GRELEASE_HBV_EC(;
    glacier_release::Number=first(@variables glacier_release),
    glacier_storage::Number=first(@variables glacier_storage),
    snow_storage::Number=first(@variables snow_storage),
    snow_liquid_storage::Number=first(@variables snow_liquid_storage),
    release_coeff::Number=first(@parameters release_coeff [description = "HBV-EC glacier-release coefficient", bounds = (0, 10), unit = "d-1"]),
    release_coeff_min::Number=first(@parameters release_coeff_min [description = "Minimum glacier-release coefficient", bounds = (0, 10), unit = "d-1"]),
    ag::Number=first(@parameters ag [description = "HBV-EC snow-cover attenuation parameter", bounds = (0, 100), unit = "mm-1"]),
    flux_name::Symbol=:grelease_hbv_ec,
)
    effective_k = release_coeff_min + (release_coeff - release_coeff_min) * exp(-ag * max(snow_storage + snow_liquid_storage, 0.0))
    @hydroflux flux_name glacier_release ~ max(0.0, effective_k * glacier_storage)
end

end
