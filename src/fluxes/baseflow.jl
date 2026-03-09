module Baseflow

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export BASE_CONSTANT,
       BASE_LINEAR,
       BASE_LINEAR_ANALYTIC,
       BASE_POWER_LAW,
       BASE_VIC,
       BASE_GR4J,
       BASE_TOPMODEL,
       BASE_THRESH_POWER,
       BASE_THRESH_STOR

"""Constant baseflow."""
function BASE_CONSTANT(;
    baseflow::Number=first(@variables baseflow),
    max_baseflow_rate::Number=first(@parameters max_baseflow_rate [description = "Maximum baseflow", bounds = (0, 100), unit = "mm/d"]),
    flux_name::Symbol=:base_constant,
)
    @hydroflux flux_name baseflow ~ max(0.0, max_baseflow_rate)
end

"""Linear storage baseflow."""
function BASE_LINEAR(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Baseflow coefficient", bounds = (0, 1), unit = "d-1"]),
    flux_name::Symbol=:base_linear,
)
    @hydroflux flux_name baseflow ~ clamp(baseflow_coeff * max(0.0, waterstorage), 0.0, max(0.0, waterstorage))
end

"""Analytical daily solution of a linear reservoir."""
function BASE_LINEAR_ANALYTIC(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Baseflow coefficient", bounds = (0, 1), unit = "d-1"]),
    flux_name::Symbol=:base_linear_analytic,
)
    @hydroflux flux_name baseflow ~ clamp(max(0.0, waterstorage) * (1 - exp(-baseflow_coeff)), 0.0, max(0.0, waterstorage))
end

"""Non-linear storage baseflow using a power law."""
function BASE_POWER_LAW(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Baseflow coefficient", bounds = (0, 100), unit = "mm^(1-n)/d"]),
    baseflow_n::Number=first(@parameters baseflow_n [description = "Baseflow exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:base_power_law,
)
    @hydroflux flux_name baseflow ~ clamp(baseflow_coeff * max(0.0, waterstorage)^baseflow_n, 0.0, max(0.0, waterstorage))
end

"""VIC baseflow with a power-law lower branch and linear upper branch."""
function BASE_VIC(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    max_baseflow_rate::Number=first(@parameters max_baseflow_rate [description = "Baseflow rate at the VIC threshold", bounds = (0, 100), unit = "mm/d"]),
    saturated_baseflow_rate::Number=first(@parameters saturated_baseflow_rate [description = "Baseflow rate at full saturation", bounds = (0, 100), unit = "mm/d"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum storage", bounds = (0, 5000), unit = "mm"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Fraction of maximum storage for the VIC breakpoint", bounds = (0, 1), unit = "-"]),
    baseflow_n::Number=first(@parameters baseflow_n [description = "VIC baseflow exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:base_vic,
)
    threshold_storage = storage_threshold * max_waterstorage
    lower_ratio = max(0.0, waterstorage) / max(threshold_storage, 1.0e-12)
    lower_branch = max_baseflow_rate * lower_ratio^baseflow_n
    upper_branch = max_baseflow_rate + (saturated_baseflow_rate - max_baseflow_rate) *
                   (max(0.0, waterstorage) - threshold_storage) /
                   max(max_waterstorage * (1 - storage_threshold), 1.0e-12)
    @hydroflux flux_name baseflow ~ clamp(
        ifelse(max(0.0, waterstorage) < threshold_storage, lower_branch, upper_branch),
        0.0,
        max(0.0, waterstorage),
    )
end

"""GR4J routing-store outflow."""
function BASE_GR4J(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    reference_waterstorage::Number=first(@parameters reference_waterstorage [description = "Reference routing storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:base_gr4j,
)
    storage_ratio = max(0.0, waterstorage) / max(reference_waterstorage, 1.0e-12)
    @hydroflux flux_name baseflow ~ clamp(
        max(0.0, waterstorage) * (1 - (1 + ((4 / 9) * storage_ratio)^4)^(-1 / 4)),
        0.0,
        max(0.0, waterstorage),
    )
end

"""TOPMODEL baseflow, where the storage variable is interpreted as a deficit."""
function BASE_TOPMODEL(;
    baseflow::Number=first(@variables baseflow),
    storage_deficit::Number=first(@variables storage_deficit),
    max_baseflow_rate::Number=first(@parameters max_baseflow_rate [description = "Baseflow rate at zero storage deficit", bounds = (0, 100), unit = "mm/d"]),
    baseflow_n::Number=first(@parameters baseflow_n [description = "TOPMODEL scaling parameter", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:base_topmodel,
)
    @hydroflux flux_name baseflow ~ max(0.0, max_baseflow_rate * exp(-max(0.0, storage_deficit) / max(baseflow_n, 1.0e-12)))
end

"""Threshold baseflow with a power-law response above the storage threshold."""
function BASE_THRESH_POWER(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Threshold baseflow coefficient", bounds = (0, 100), unit = "mm^(1-n)/d"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Storage threshold above which baseflow begins", bounds = (0, 5000), unit = "mm"]),
    baseflow_n::Number=first(@parameters baseflow_n [description = "Threshold baseflow exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:base_thresh_power,
)
    @hydroflux flux_name baseflow ~ clamp(
        baseflow_coeff * max(max(0.0, waterstorage) - storage_threshold, 0.0)^baseflow_n,
        0.0,
        max(0.0, waterstorage),
    )
end

"""Threshold baseflow with a linear response above the storage threshold."""
function BASE_THRESH_STOR(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Baseflow coefficient", bounds = (0, 100), unit = "d-1"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Storage threshold above which baseflow begins", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:base_thresh_stor,
)
    @hydroflux flux_name baseflow ~ clamp(baseflow_coeff * max(max(0.0, waterstorage) - storage_threshold, 0.0), 0.0, max(0.0, waterstorage))
end

end
