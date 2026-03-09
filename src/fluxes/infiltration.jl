module Infiltration

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export INF_RATIONAL,
       INF_SCS,
       INF_GREEN_AMPT,
       INF_GA_SIMPLE,
       INF_VIC,
       INF_VIC_ARNO,
       INF_HBV,
       INF_PRMS,
       INF_UBC,
       INF_GR4J,
       INF_HMETS,
       INF_AWBM,
       INF_PDM

_pdm_cstar(storage::Number, max_storage::Number, pdm_b::Number) =
    max_storage * (pdm_b + 1) * (1 - (1 - clamp(max(0.0, storage) / max(max_storage, 1.0e-12), 0.0, 1.0))^(1 / (pdm_b + 1)))

_pdm_storage_from_cstar(cstar::Number, max_storage::Number, pdm_b::Number) =
    max_storage * (1 - (1 - clamp(cstar / max(max_storage * (pdm_b + 1), 1.0e-12), 0.0, 1.0))^(pdm_b + 1))

"""Rational runoff partition infiltration."""
function INF_RATIONAL(;
    infiltration::Number=first(@variables infiltration),
    rainfall::Number=first(@variables rainfall),
    runoff_coeff::Number=first(@parameters runoff_coeff [description = "Runoff coefficient", bounds = (0, 1), unit = "-"]),
    flux_name::Symbol=:inf_rational,
)
    @hydroflux flux_name infiltration ~ max(0.0, rainfall) * (1 - clamp(runoff_coeff, 0.0, 1.0))
end

"""SCS-CN infiltration."""
function INF_SCS(;
    infiltration::Number=first(@variables infiltration),
    rainfall::Number=first(@variables rainfall),
    curve_number::Number=first(@parameters curve_number [description = "SCS curve number", bounds = (1, 100), unit = "-"]),
    flux_name::Symbol=:inf_scs,
)
    s = max(25400 / max(curve_number, 1.0e-12) - 254, 0.0)
    runoff = ifelse(rainfall <= 0.2 * s, 0.0, (max(0.0, rainfall - 0.2 * s)^2) / max(rainfall + 0.8 * s, 1.0e-12))
    @hydroflux flux_name infiltration ~ clamp(max(0.0, rainfall) - runoff, 0.0, max(0.0, rainfall))
end

"""Green-Ampt infiltration."""
function INF_GREEN_AMPT(;
    infiltration::Number=first(@variables infiltration),
    rainfall::Number=first(@variables rainfall),
    sat_hydraulic_cond::Number=first(@parameters sat_hydraulic_cond [description = "Saturated hydraulic conductivity", bounds = (0, 5000), unit = "mm/d"]),
    wetting_front_suction::Number=first(@parameters wetting_front_suction [description = "Wetting-front suction head", bounds = (0, 5000), unit = "mm"]),
    max_porosity::Number=first(@parameters max_porosity [description = "Maximum porosity", bounds = (0, 1), unit = "-"]),
    soil_porosity::Number=first(@variables soil_porosity),
    cumulative_infiltration::Number=first(@variables cumulative_infiltration),
    flux_name::Symbol=:inf_green_ampt,
)
    capacity = sat_hydraulic_cond * (1 + wetting_front_suction * max(max_porosity - soil_porosity, 0.0) / max(cumulative_infiltration, 1.0e-12))
    @hydroflux flux_name infiltration ~ min(max(0.0, rainfall), max(0.0, capacity))
end

INF_GA_SIMPLE(; kwargs...) = INF_GREEN_AMPT(; kwargs..., flux_name=get(kwargs, :flux_name, :inf_ga_simple))

"""VIC infiltration."""
function INF_VIC(;
    infiltration::Number=first(@variables infiltration),
    rainfall::Number=first(@variables rainfall),
    soil_storage::Number=first(@variables soil_storage),
    max_soil_storage::Number=first(@parameters max_soil_storage [description = "Maximum soil storage", bounds = (0, 5000), unit = "mm"]),
    vic_gamma::Number=first(@parameters vic_gamma [description = "VIC shape parameter", bounds = (0.01, 20), unit = "-"]),
    vic_alpha::Number=first(@parameters vic_alpha [description = "VIC alpha parameter", bounds = (0.01, 20), unit = "-"]),
    zmax::Number=first(@parameters zmax [description = "Maximum infiltration depth", bounds = (0.01, 5000), unit = "mm"]),
    zmin::Number=first(@parameters zmin [description = "Minimum infiltration depth", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:inf_vic,
)
    k1 = ((max(zmax - zmin, 1.0e-12)) * vic_alpha * vic_gamma)^(-vic_gamma)
    storage_ratio = clamp(max(0.0, soil_storage) / max(max_soil_storage, 1.0e-12), 0.0, 1.0)
    inner = max(vic_gamma * vic_alpha * zmax + zmin - storage_ratio, 0.0)
    @hydroflux flux_name infiltration ~ clamp(max(0.0, rainfall) * k1 * inner^vic_gamma, 0.0, max(0.0, rainfall))
end

"""VIC/ARNO infiltration."""
function INF_VIC_ARNO(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    exp_coeff::Number=first(@parameters exp_coeff [description = "ARNO infiltration exponent", bounds = (0, 10), unit = "-"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:inf_vic_arno,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name infiltration ~ max(0.0, prcpeff) * (1 - (1 - storage_ratio)^exp_coeff)
end

"""HBV infiltration."""
function INF_HBV(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    hbv_beta::Number=first(@parameters hbv_beta [description = "HBV non-linearity parameter", bounds = (0, 10), unit = "-"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:inf_hbv,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name infiltration ~ max(0.0, prcpeff) * (1 - storage_ratio^hbv_beta)
end

"""PRMS infiltration as rainfall reduced by the saturated contributing area."""
function INF_PRMS(;
    infiltration::Number=first(@variables infiltration),
    rainfall::Number=first(@variables rainfall),
    soil_storage::Number=first(@variables soil_storage),
    tension_storage::Number=first(@parameters tension_storage [description = "Tension storage", bounds = (0, 5000), unit = "mm"]),
    fsat_max::Number=first(@parameters fsat_max [description = "Maximum saturated contributing area", bounds = (0, 1), unit = "-"]),
    flux_name::Symbol=:inf_prms,
)
    saturation_fraction = fsat_max * min(max(0.0, soil_storage) / max(tension_storage, 1.0e-12), 1.0)
    @hydroflux flux_name infiltration ~ max(0.0, rainfall) * (1 - clamp(saturation_fraction, 0.0, 1.0))
end

"""UBC watershed model infiltration partition."""
function INF_UBC(;
    infiltration::Number=first(@variables infiltration),
    percolation::Number=first(@variables percolation),
    interflow::Number=first(@variables interflow),
    runoff::Number=first(@variables runoff),
    rainfall::Number=first(@variables rainfall),
    soil_storage::Number=first(@variables soil_storage),
    max_soil_storage::Number=first(@parameters max_soil_storage [description = "Maximum soil storage", bounds = (0, 5000), unit = "mm"]),
    ponded_water::Number=first(@variables ponded_water),
    max_perc_rate::Number=first(@parameters max_perc_rate [description = "Maximum percolation rate", bounds = (0, 5000), unit = "mm/d"]),
    imp_fraction::Number=first(@parameters imp_fraction [description = "Impervious fraction", bounds = (0, 1), unit = "-"]),
    poagfn::Number=first(@parameters poagfn [description = "UBC infiltration-shape parameter", bounds = (1.0e-6, 5000), unit = "mm"]),
    voflax::Number=first(@parameters voflax [description = "UBC ponding parameter", bounds = (1.0e-6, 5000), unit = "mm"]),
    infil_flux_name::Symbol=:inf_ubc,
    perc_flux_name::Symbol=:inf_ubc_perc,
    interflow_flux_name::Symbol=:inf_ubc_interflow,
    runoff_flux_name::Symbol=:inf_ubc_runoff,
)
    b1 = imp_fraction * 10^(-(max_soil_storage - soil_storage) / max(poagfn, 1.0e-12))
    ff = 1 + log(max(ponded_water, 1.0e-12) / max(voflax, 1.0e-12)) / log(max(voflax, 1.0e-12) / 1800)
    b2 = clamp(b1 + (1 - b1) * ff, 0.0, 1.0)
    inf_val = max(0.0, rainfall) * (1 - b2)
    excess = max(max(0.0, rainfall) - inf_val, 0.0)
    perc_val = min(max_perc_rate, excess) * (1 - b2)
    interflow_val = max(excess - perc_val, 0.0) * (1 - b2)
    runoff_val = b2 * max(0.0, rainfall)
    infil_flux = @hydroflux infil_flux_name infiltration ~ inf_val
    perc_flux = @hydroflux perc_flux_name percolation ~ perc_val
    inter_flux = @hydroflux interflow_flux_name interflow ~ interflow_val
    runoff_flux = @hydroflux runoff_flux_name runoff ~ runoff_val
    return (infil_flux, perc_flux, inter_flux, runoff_flux)
end

"""GR4J production-store infiltration."""
function INF_GR4J(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Production store capacity", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:inf_gr4j,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    tanh_term = tanh(max(0.0, prcpeff) / max(max_waterstorage, 1.0e-12))
    @hydroflux flux_name infiltration ~ clamp(
        max_waterstorage * (1 - storage_ratio^2) * tanh_term / (1 + storage_ratio * tanh_term),
        0.0,
        max(0.0, prcpeff),
    )
end

"""HMETS infiltration."""
function INF_HMETS(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    runoff_coeff::Number=first(@parameters runoff_coeff [description = "Runoff coefficient", bounds = (0, 1), unit = "-"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum storage", bounds = (0, 5000), unit = "mm"]),
    flux_name::Symbol=:inf_hmets,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name infiltration ~ max(0.0, prcpeff) * clamp(1 - runoff_coeff * storage_ratio, 0.0, 1.0)
end

"""AWBM infiltration derived from three partial-area stores."""
function INF_AWBM(;
    infiltration::Number=first(@variables infiltration),
    excess::Number=first(@variables excess),
    runoff::Number=first(@variables runoff),
    to_groundwater::Number=first(@variables to_groundwater),
    rainfall::Number=first(@variables rainfall),
    store1::Number=first(@variables store1),
    store2::Number=first(@variables store2),
    store3::Number=first(@variables store3),
    area1::Number=first(@parameters area1 [description = "Area fraction 1", bounds = (0, 1), unit = "-"]),
    area2::Number=first(@parameters area2 [description = "Area fraction 2", bounds = (0, 1), unit = "-"]),
    area3::Number=first(@parameters area3 [description = "Area fraction 3", bounds = (0, 1), unit = "-"]),
    max_store1::Number=first(@parameters max_store1 [description = "Maximum store 1 capacity", bounds = (0, 5000), unit = "mm"]),
    max_store2::Number=first(@parameters max_store2 [description = "Maximum store 2 capacity", bounds = (0, 5000), unit = "mm"]),
    max_store3::Number=first(@parameters max_store3 [description = "Maximum store 3 capacity", bounds = (0, 5000), unit = "mm"]),
    bfi::Number=first(@parameters bfi [description = "Baseflow index", bounds = (0, 1), unit = "-"]),
    infil_flux_name::Symbol=:inf_awbm,
    excess_flux_name::Symbol=:inf_awbm_excess,
    runoff_flux_name::Symbol=:inf_awbm_runoff,
    gw_flux_name::Symbol=:inf_awbm_to_gw,
)
    excess1 = max(area1 * rainfall - max(area1 * max_store1 - store1, 0.0), 0.0)
    excess2 = max(area2 * rainfall - max(area2 * max_store2 - store2, 0.0), 0.0)
    excess3 = max(area3 * rainfall - max(area3 * max_store3 - store3, 0.0), 0.0)
    excess_val = excess1 + excess2 + excess3
    infil_val = max(max(0.0, rainfall) - excess_val, 0.0)
    runoff_val = (1 - clamp(bfi, 0.0, 1.0)) * excess_val
    gw_val = clamp(bfi, 0.0, 1.0) * excess_val
    infil_flux = @hydroflux infil_flux_name infiltration ~ infil_val
    excess_flux = @hydroflux excess_flux_name excess ~ excess_val
    runoff_flux = @hydroflux runoff_flux_name runoff ~ runoff_val
    gw_flux = @hydroflux gw_flux_name to_groundwater ~ gw_val
    return (infil_flux, excess_flux, runoff_flux, gw_flux)
end

"""PDM infiltration inferred from the standard probability-distributed storage relation."""
function INF_PDM(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Mean storage capacity", bounds = (0, 5000), unit = "mm"]),
    pdm_b::Number=first(@parameters pdm_b [description = "PDM shape parameter", bounds = (0, 10), unit = "-"]),
    flux_name::Symbol=:inf_pdm,
)
    current_cstar = _pdm_cstar(waterstorage, max_waterstorage, pdm_b)
    cmax = max_waterstorage * (pdm_b + 1)
    next_cstar = min(cmax, current_cstar + max(0.0, prcpeff))
    next_storage = _pdm_storage_from_cstar(next_cstar, max_waterstorage, pdm_b)
    @hydroflux flux_name infiltration ~ clamp(next_storage - max(0.0, waterstorage), 0.0, max(0.0, prcpeff))
end

end
