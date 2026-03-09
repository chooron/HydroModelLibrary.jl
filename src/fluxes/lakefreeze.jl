module LakeFreeze

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export LFREEZE_BASIC

function LFREEZE_BASIC(;
    lake_freezing::Number=first(@variables lake_freezing),
    snow_water_equivalent::Number=first(@variables snow_water_equivalent),
    potential_melt::Number=first(@variables potential_melt),
    freeze_scale::Number=first(@parameters freeze_scale [description = "Snow insulation scale", bounds = (1.0e-6, 5000), unit = "mm"]),
    flux_name::Symbol=:lfreeze_basic,
)
    insulation = 1 - min(max(0.0, snow_water_equivalent) / max(freeze_scale, 1.0e-12), 1.0)
    @hydroflux flux_name lake_freezing ~ -insulation * potential_melt
end

end
