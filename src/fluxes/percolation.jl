module Percolation

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export PERC_CONSTANT,
       PERC_LINEAR,
       PERC_POWER_LAW,
       PERC_PRMS,
       PERC_SACRAMENTO,
       PERC_GAWSER,
       PERC_GAWSER_CONSTRAIN,
       PERC_GR4JEXCH,
       PERC_GR4JEXCH2

"""Constant percolation limited by available storage."""
function PERC_CONSTANT(;
    percolation::Number=first(@variables percolation),
    waterstorage::Number=first(@variables waterstorage),
    max_perc_rate::Number=first(@parameters max_perc_rate [description = "Maximum percolation rate", bounds = (0, 5000), unit = "mm/d"]),
    flux_name::Symbol=:perc_constant,
)
    @hydroflux flux_name percolation ~ min(max(0.0, waterstorage), max(0.0, max_perc_rate))
end

"""Linear percolation."""
function PERC_LINEAR(;
    percolation::Number=first(@variables percolation),
    waterstorage::Number=first(@variables waterstorage),
    perc_coeff::Number=first(@parameters perc_coeff [description = "Linear percolation coefficient", bounds = (0, 10), unit = "d-1"]),
    flux_name::Symbol=:perc_linear,
)
    @hydroflux flux_name percolation ~ clamp(perc_coeff * max(0.0, waterstorage), 0.0, max(0.0, waterstorage))
end

"""Power-law percolation."""
function PERC_POWER_LAW(;
    percolation::Number=first(@variables percolation),
    waterstorage::Number=first(@variables waterstorage),
    perc_coeff::Number=first(@parameters perc_coeff [description = "Power-law percolation coefficient", bounds = (0, 100), unit = "mm^(1-n)/d"]),
    perc_n::Number=first(@parameters perc_n [description = "Power-law percolation exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:perc_power_law,
)
    @hydroflux flux_name percolation ~ clamp(perc_coeff * max(0.0, waterstorage)^perc_n, 0.0, max(0.0, waterstorage))
end

"""PRMS-style gravity drainage / recharge."""
function PERC_PRMS(;
    percolation::Number=first(@variables percolation),
    waterstorage::Number=first(@variables waterstorage),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum storage", bounds = (0, 5000), unit = "mm"]),
    sat_hydraulic_cond::Number=first(@parameters sat_hydraulic_cond [description = "Saturated hydraulic conductivity", bounds = (0, 5000), unit = "mm/d"]),
    perc_n::Number=first(@parameters perc_n [description = "PRMS recharge exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:perc_prms,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name percolation ~ clamp(sat_hydraulic_cond * storage_ratio^perc_n, 0.0, max(0.0, waterstorage))
end

"""Sacramento lower-zone percolation demand formulation."""
function PERC_SACRAMENTO(;
    percolation::Number=first(@variables percolation),
    upper_waterstorage::Number=first(@variables upper_waterstorage),
    upper_waterstorage_max::Number=first(@parameters upper_waterstorage_max [description = "Upper-zone free-water capacity", bounds = (0, 5000), unit = "mm"]),
    lower_zone_deficit::Number=first(@variables lower_zone_deficit),
    pbase::Number=first(@parameters pbase [description = "Base percolation demand", bounds = (0, 5000), unit = "mm/d"]),
    zperc::Number=first(@parameters zperc [description = "Sacramento deficit scaling coefficient", bounds = (0, 100), unit = "-"]),
    rexp::Number=first(@parameters rexp [description = "Sacramento deficit exponent", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:perc_sacramento,
)
    demand = pbase * (1 + zperc * max(0.0, lower_zone_deficit)^(1 + rexp))
    @hydroflux flux_name percolation ~ clamp(demand * max(0.0, upper_waterstorage) / max(upper_waterstorage_max, 1.0e-12), 0.0, max(0.0, upper_waterstorage))
end

"""GAWSER threshold percolation."""
function PERC_GAWSER(;
    percolation::Number=first(@variables percolation),
    waterstorage::Number=first(@variables waterstorage),
    perc_coeff::Number=first(@parameters perc_coeff [description = "GAWSER percolation coefficient", bounds = (0, 10), unit = "d-1"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Storage threshold", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:perc_gawser,
)
    @hydroflux flux_name percolation ~ clamp(perc_coeff * max(max(0.0, waterstorage) - storage_threshold, 0.0), 0.0, max(0.0, waterstorage))
end

"""GAWSER threshold percolation constrained by an explicit available storage."""
function PERC_GAWSER_CONSTRAIN(;
    percolation::Number=first(@variables percolation),
    waterstorage::Number=first(@variables waterstorage),
    available_storage::Number=first(@variables available_storage),
    perc_coeff::Number=first(@parameters perc_coeff [description = "GAWSER percolation coefficient", bounds = (0, 10), unit = "d-1"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "Storage threshold", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:perc_gawser_constrain,
)
    raw_perc = perc_coeff * max(max(0.0, waterstorage) - storage_threshold, 0.0)
    @hydroflux flux_name percolation ~ clamp(raw_perc, 0.0, max(0.0, available_storage))
end

"""GR4J groundwater exchange term expressed on the routing-store state."""
function PERC_GR4JEXCH(;
    percolation::Number=first(@variables percolation),
    waterstorage::Number=first(@variables waterstorage),
    x2::Number=first(@parameters x2 [description = "GR4J exchange coefficient", bounds = (-5000, 5000), unit = "mm/d"]),
    x3::Number=first(@parameters x3 [description = "GR4J routing-store reference capacity", bounds = (1.0e-6, 5000), unit = "mm"]),
    flux_name::Symbol=:perc_gr4jexch,
)
    storage_ratio = min(max(max(0.0, waterstorage) / max(x3, 1.0e-12), 0.0), 1.0)
    @hydroflux flux_name percolation ~ -x2 * storage_ratio^3.5
end

"""Alternative GR4J exchange alias used by Raven."""
function PERC_GR4JEXCH2(; kwargs...)
    params = (; kwargs...)
    return PERC_GR4JEXCH(; params..., flux_name=get(params, :flux_name, :perc_gr4jexch2))
end

end
