module GlacierInfiltration

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export GLACIERINF_ALL,
       GLACIERINF_LINEAR,
       GLACIERINF_GSMSOCONT

"""Direct transfer of all available liquid water into a glacier routing store."""
function GLACIERINF_ALL(;
    glacier_infiltration::Number=first(@variables glacier_infiltration),
    liquid_water::Number=first(@variables liquid_water),
    flux_name::Symbol=:glacierinf_all,
)
    @hydroflux flux_name glacier_infiltration ~ max(0.0, liquid_water)
end

"""Linear glacier infiltration as a fraction of available liquid water."""
function GLACIERINF_LINEAR(;
    glacier_infiltration::Number=first(@variables glacier_infiltration),
    liquid_water::Number=first(@variables liquid_water),
    infil_fraction::Number=first(@parameters infil_fraction [description = "Fraction of available liquid water infiltrating into glacier storage", bounds = (0, 1), unit = "-"]),
    flux_name::Symbol=:glacierinf_linear,
)
    @hydroflux flux_name glacier_infiltration ~ clamp(infil_fraction, 0.0, 1.0) * max(0.0, liquid_water)
end

"""GSMSOCONT-style rain-on-ice input that is suppressed when a snowpack covers the glacier."""
function GLACIERINF_GSMSOCONT(;
    glacier_infiltration::Number=first(@variables glacier_infiltration),
    rain_on_ice::Number=first(@variables rain_on_ice),
    snow_storage::Number=first(@variables snow_storage),
    flux_name::Symbol=:glacierinf_gsmsocont,
)
    @hydroflux flux_name glacier_infiltration ~ ifelse(snow_storage > 0.0, 0.0, max(0.0, rain_on_ice))
end

end
