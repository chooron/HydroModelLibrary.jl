module Abstraction

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export ABST_FILL,
       ABST_MAX,
       ABST_RATIO,
       ABST_DERIVED,
       ABST_SCS,
       ABST_PERCENTAGE,
       ABST_PDMROF,
       ABST_UWFS

"""Fill abstraction storage up to a target maximum."""
function ABST_FILL(;
    abstraction::Number=first(@variables abstraction),
    precipitation::Number=first(@variables precipitation),
    waterstorage::Number=first(@variables waterstorage),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum abstraction storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:abst_fill,
)
    available_capacity = max(max_waterstorage - max(0.0, waterstorage), 0.0)
    @hydroflux flux_name abstraction ~ min(max(0.0, precipitation), available_capacity)
end

"""Maximum abstraction rate."""
function ABST_MAX(;
    abstraction::Number=first(@variables abstraction),
    precipitation::Number=first(@variables precipitation),
    max_abstraction::Number=first(@parameters max_abstraction [description = "Maximum abstraction rate", bounds = (0, 5000), unit = "mm/d"]),
    flux_name::Symbol=:abst_max,
)
    @hydroflux flux_name abstraction ~ min(max(0.0, precipitation), max(0.0, max_abstraction))
end

"""Fixed-ratio abstraction."""
function ABST_RATIO(;
    abstraction::Number=first(@variables abstraction),
    precipitation::Number=first(@variables precipitation),
    abstraction_ratio::Number=first(@parameters abstraction_ratio [description = "Abstraction ratio", bounds = (0, 1), unit = "-"]),
    flux_name::Symbol=:abst_ratio,
)
    @hydroflux flux_name abstraction ~ max(0.0, precipitation) * clamp(abstraction_ratio, 0.0, 1.0)
end

"""Residual abstraction after other partitioning fluxes are applied."""
function ABST_DERIVED(;
    abstraction::Number=first(@variables abstraction),
    precipitation::Number=first(@variables precipitation),
    outgoing_fluxes::Number=first(@variables outgoing_fluxes),
    flux_name::Symbol=:abst_derived,
)
    @hydroflux flux_name abstraction ~ max(max(0.0, precipitation) - max(0.0, outgoing_fluxes), 0.0)
end

"""SCS initial abstraction based on curve number."""
function ABST_SCS(;
    abstraction::Number=first(@variables abstraction),
    ponded_water::Number=first(@variables ponded_water),
    curve_number::Number=first(@parameters curve_number [description = "SCS curve number", bounds = (1, 100), unit = "-"]),
    scs_fraction::Number=first(@parameters scs_fraction [description = "Fraction of SCS retention used for abstraction", bounds = (0, 1), unit = "-"]),
    dt::Number=first(@variables dt [description = "Time step in days", bounds = (1.0e-6, 365)]),
    flux_name::Symbol=:abst_scs,
)
    retention = max(scs_fraction * 25.4 * (1000 / max(curve_number, 1.0e-12) - 10), 0.0)
    @hydroflux flux_name abstraction ~ max(retention, max(0.0, ponded_water)) / max(dt, 1.0e-12)
end

"""Percentage abstraction of ponded water."""
function ABST_PERCENTAGE(;
    abstraction::Number=first(@variables abstraction),
    ponded_water::Number=first(@variables ponded_water),
    abstraction_ratio::Number=first(@parameters abstraction_ratio [description = "Percentage abstraction", bounds = (0, 1), unit = "-"]),
    flux_name::Symbol=:abst_percentage,
)
    @hydroflux flux_name abstraction ~ clamp(abstraction_ratio, 0.0, 1.0) * max(0.0, ponded_water)
end

"""PDMROF abstraction using the probability-distributed storage relation."""
function ABST_PDMROF(;
    abstraction::Number=first(@variables abstraction),
    ponded_water::Number=first(@variables ponded_water),
    depression_storage::Number=first(@variables depression_storage),
    max_storage::Number=first(@parameters max_storage [description = "Mean depression storage", bounds = (0, 5000), unit = "mm"]),
    pdm_b::Number=first(@parameters pdm_b [description = "PDM shape parameter", bounds = (0, 20), unit = "-"]),
    dt::Number=first(@variables dt [description = "Time step in days", bounds = (1.0e-6, 365)]),
    flux_name::Symbol=:abst_pdmrof,
)
    cmax = max_storage * (pdm_b + 1)
    cstar = cmax * (1 - (1 - clamp(max(0.0, depression_storage) / max(max_storage, 1.0e-12), 0.0, 1.0))^(1 / (pdm_b + 1)))
    term1 = (1 - clamp(cstar / max(cmax, 1.0e-12), 0.0, 1.0))^(pdm_b + 1)
    term2 = (1 - clamp((cstar + max(0.0, ponded_water)) / max(cmax, 1.0e-12), 0.0, 1.0))^(pdm_b + 1)
    @hydroflux flux_name abstraction ~ max_storage * max(term1 - term2, 0.0) / max(dt, 1.0e-12)
end

"""UWFS abstraction as filling of a running storage deficit."""
function ABST_UWFS(;
    abstraction::Number=first(@variables abstraction),
    precipitation::Number=first(@variables precipitation),
    deficit_storage::Number=first(@variables deficit_storage),
    beta_ave::Number=first(@parameters beta_ave [description = "Average filling coefficient", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:abst_uwfs,
)
    @hydroflux flux_name abstraction ~ min(beta_ave * max(0.0, precipitation), max(0.0, deficit_storage))
end

end
