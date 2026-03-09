module SoilBalance

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export SOILBAL_BUCKET,
       SOILBAL_PRMS,
       SOILBAL_PENMAN

"""Generic single-store soil-water balance bucket."""
function SOILBAL_BUCKET(;
    soil_storage::Number=first(@variables soil_storage),
    infiltration::Number=first(@variables infiltration),
    soil_evaporation::Number=first(@variables soil_evaporation),
    percolation::Number=first(@variables percolation),
    overflow::Number=first(@variables overflow),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        dfluxes = begin
            @stateflux soil_storage ~ infiltration - soil_evaporation - percolation - overflow
        end
    end
end

"""PRMS-style two-layer soil-balance bucket with recharge and lower soil stores."""
function SOILBAL_PRMS(;
    recharge_storage::Number=first(@variables recharge_storage),
    lower_storage::Number=first(@variables lower_storage),
    infiltration::Number=first(@variables infiltration),
    upper_evaporation::Number=first(@variables upper_evaporation),
    percolation::Number=first(@variables percolation),
    lower_evaporation::Number=first(@variables lower_evaporation),
    excess::Number=first(@variables excess),
    lower_storage_max::Number=first(@parameters lower_storage_max [description = "Maximum lower soil storage", bounds = (0, 5000), unit = "mm"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux excess ~ max(lower_storage + percolation - lower_storage_max, 0.0)
        end
        dfluxes = begin
            @stateflux recharge_storage ~ infiltration - upper_evaporation - percolation
            @stateflux lower_storage ~ percolation - lower_evaporation - excess
        end
    end
end

"""Penman-style root-zone storage plus deficit store balance."""
function SOILBAL_PENMAN(;
    rootzone_storage::Number=first(@variables rootzone_storage),
    deficit_storage::Number=first(@variables deficit_storage),
    precipitation::Number=first(@variables precipitation),
    direct_evaporation::Number=first(@variables direct_evaporation),
    saturation_excess::Number=first(@variables saturation_excess),
    transfer_to_deficit::Number=first(@variables transfer_to_deficit),
    transpiration::Number=first(@variables transpiration),
    deficit_release::Number=first(@variables deficit_release),
    smax::Number=first(@parameters smax [description = "Maximum root-zone storage", bounds = (0, 5000), unit = "mm"]),
    split_fraction::Number=first(@parameters split_fraction [description = "Fraction of saturation excess routed directly", bounds = (0, 1), unit = "-"]),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux saturation_excess ~ max(rootzone_storage + precipitation - smax, 0.0)
            @hydroflux transfer_to_deficit ~ (1 - split_fraction) * saturation_excess
            @hydroflux deficit_release ~ min(max(0.0, deficit_storage), max(0.0, transfer_to_deficit))
        end
        dfluxes = begin
            @stateflux rootzone_storage ~ precipitation - direct_evaporation - saturation_excess
            @stateflux deficit_storage ~ transfer_to_deficit - transpiration - deficit_release
        end
    end
end

end
