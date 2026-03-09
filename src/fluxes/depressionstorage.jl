module DepressionStorage

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export DEPSTOR_BUCKET,
       DEPSTOR_PRMS,
       DFLOW_THRESHPOW,
       DFLOW_LINEAR,
       DFLOW_WEIR

"""Generic depression or wetland storage bucket."""
function DEPSTOR_BUCKET(;
    depression_storage::Number=first(@variables depression_storage),
    trapped_inflow::Number=first(@variables trapped_inflow),
    evaporation::Number=first(@variables evaporation),
    seepage::Number=first(@variables seepage),
    overflow::Number=first(@variables overflow),
    name::Union{Nothing, Symbol}=nothing,
)
    return HydroModels.@hydrobucket name begin
        dfluxes = begin
            @stateflux depression_storage ~ trapped_inflow - evaporation - seepage - overflow
        end
    end
end

"""PRMS/impervious-style depression storage overflow once capacity is exceeded."""
function DEPSTOR_PRMS(;
    overflow::Number=first(@variables overflow),
    depression_storage::Number=first(@variables depression_storage),
    incoming_water::Number=first(@variables incoming_water),
    max_depression_storage::Number=first(@parameters max_depression_storage [description = "Maximum depression storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:depstor_prms,
)
    @hydroflux flux_name overflow ~ min(max(0.0, incoming_water), max(max(0.0, depression_storage) - max_depression_storage, 0.0))
end

"""Threshold-power depression or wetland overflow."""
function DFLOW_THRESHPOW(;
    depression_flow::Number=first(@variables depression_flow),
    depression_storage::Number=first(@variables depression_storage),
    max_flow::Number=first(@parameters max_flow [description = "Maximum depression flow", bounds = (0, 5000), unit = "mm/d"]),
    max_storage::Number=first(@parameters max_storage [description = "Maximum depression storage", bounds = (0, 5000), unit = "mm"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Depression flow threshold ratio", bounds = (0, 1), unit = "-"]),
    flow_n::Number=first(@parameters flow_n [description = "Depression-flow exponent", bounds = (0, 20), unit = "-"]),
    flux_name::Symbol=:dflow_threshpow,
)
    storage_ratio = clamp(max(0.0, depression_storage) / max(max_storage, 1.0e-12), 0.0, 1.0)
    response = clamp((storage_ratio - storage_threshold) / max(1 - storage_threshold, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name depression_flow ~ max_flow * response^flow_n
end

"""Linear depression or wetland overflow above a threshold."""
function DFLOW_LINEAR(;
    depression_flow::Number=first(@variables depression_flow),
    depression_storage::Number=first(@variables depression_storage),
    drain_coeff::Number=first(@parameters drain_coeff [description = "Depression drainage coefficient", bounds = (0, 100), unit = "d-1"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Depression storage threshold", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:dflow_linear,
)
    @hydroflux flux_name depression_flow ~ max(0.0, drain_coeff * max(depression_storage - storage_threshold, 0.0))
end

"""Weir-flow depression overflow."""
function DFLOW_WEIR(;
    depression_flow::Number=first(@variables depression_flow),
    depression_storage::Number=first(@variables depression_storage),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Weir threshold storage", bounds = (0, 5000), unit = "mm"]),
    depression_ratio::Number=first(@parameters depression_ratio [description = "Depression geometry ratio", bounds = (0, 1000), unit = "-"]),
    area::Number=first(@parameters area [description = "Depression area", bounds = (1.0e-6, 1.0e12), unit = "m2"]),
    flux_name::Symbol=:dflow_weir,
)
    g = 9.80665
    coeff = depression_ratio * sqrt(max(area, 1.0e-12)) / max(area, 1.0e-12)
    head = max(depression_storage - storage_threshold, 0.0)
    @hydroflux flux_name depression_flow ~ 0.666 * coeff * sqrt(2 * g) * head^1.5
end

end
