module SoilEvap

using ..HydroModels
using ..HydroModels: @parameters, @variables, Number

export SOILEVAP_ALL,
       SOILEVAP_LINEAR,
       SOILEVAP_ROOT,
       SOILEVAP_SEQUEN,
       SOILEVAP_TOPMODEL,
       SOILEVAP_VIC,
       SOILEVAP_HBV,
       SOILEVAP_HBV_ORESUND,
       SOILEVAP_PRMS,
       SOILEVAP_SACSMA,
       SOILEVAP_GR4J,
       SOILEVAP_UBC,
       SOILEVAP_PDM,
       SOILEVAP_HYPR,
       SOILEVAP_CHU,
       SOILEVAP_AWBM,
       SOILEVAP_HYMOD2

function SOILEVAP_ALL(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    potential_evaporation::Number=first(@variables potential_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_all,
)
    @hydroflux flux_name soil_evaporation ~ clamp(potential_evaporation * pet_corr, 0.0, max(0.0, waterstorage))
end

function SOILEVAP_LINEAR(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    evap_coeff::Number=first(@parameters evap_coeff [description = "Linear evaporation coefficient", bounds = (0, 100), unit = "d-1"]),
    flux_name::Symbol=:soilevap_linear,
)
    @hydroflux flux_name soil_evaporation ~ min(max(0.0, potential_evaporation), max(0.0, evap_coeff * waterstorage))
end

function SOILEVAP_ROOT(;
    upper_evaporation::Number=first(@variables upper_evaporation),
    lower_evaporation::Number=first(@variables lower_evaporation),
    upper_storage::Number=first(@variables upper_storage),
    lower_storage::Number=first(@variables lower_storage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    upper_tension::Number=first(@parameters upper_tension [description = "Upper tension storage", bounds = (0, 5000), unit = "mm"]),
    lower_tension::Number=first(@parameters lower_tension [description = "Lower tension storage", bounds = (0, 5000), unit = "mm"]),
    upper_root_fraction::Number=first(@parameters upper_root_fraction [description = "Upper-root fraction", bounds = (0, 1), unit = "-"]),
    lower_root_fraction::Number=first(@parameters lower_root_fraction [description = "Lower-root fraction", bounds = (0, 1), unit = "-"]),
    upper_flux_name::Symbol=:soilevap_root_upper,
    lower_flux_name::Symbol=:soilevap_root_lower,
)
    upper_ratio = min(max(0.0, upper_storage) / max(upper_tension, 1.0e-12), 1.0)
    lower_ratio = min(max(0.0, lower_storage) / max(lower_tension, 1.0e-12), 1.0)
    upper_flux = @hydroflux upper_flux_name upper_evaporation ~ max(0.0, potential_evaporation) * upper_root_fraction * upper_ratio
    lower_flux = @hydroflux lower_flux_name lower_evaporation ~ max(0.0, potential_evaporation) * lower_root_fraction * lower_ratio
    return (upper_flux, lower_flux)
end

function SOILEVAP_SEQUEN(;
    upper_evaporation::Number=first(@variables upper_evaporation),
    lower_evaporation::Number=first(@variables lower_evaporation),
    upper_storage::Number=first(@variables upper_storage),
    lower_storage::Number=first(@variables lower_storage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    upper_tension::Number=first(@parameters upper_tension [description = "Upper tension storage", bounds = (0, 5000), unit = "mm"]),
    lower_tension::Number=first(@parameters lower_tension [description = "Lower tension storage", bounds = (0, 5000), unit = "mm"]),
    upper_flux_name::Symbol=:soilevap_sequen_upper,
    lower_flux_name::Symbol=:soilevap_sequen_lower,
)
    upper_val = max(0.0, potential_evaporation) * min(max(0.0, upper_storage) / max(upper_tension, 1.0e-12), 1.0)
    lower_val = max(0.0, potential_evaporation - upper_val) * min(max(0.0, lower_storage) / max(lower_tension, 1.0e-12), 1.0)
    upper_flux = @hydroflux upper_flux_name upper_evaporation ~ upper_val
    lower_flux = @hydroflux lower_flux_name lower_evaporation ~ lower_val
    return (upper_flux, lower_flux)
end

function SOILEVAP_TOPMODEL(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    storage_tension::Number=first(@parameters storage_tension [description = "Tension storage", bounds = (0, 5000), unit = "mm"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_topmodel,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(storage_tension, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name soil_evaporation ~ clamp(potential_evaporation * pet_corr * storage_ratio, 0.0, max(0.0, waterstorage))
end

function SOILEVAP_VIC(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum soil storage", bounds = (0, 5000), unit = "mm"]),
    evap_gamma::Number=first(@parameters evap_gamma [description = "VIC evaporation exponent", bounds = (0.01, 20), unit = "-"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_vic,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    evap_factor = 1 - (1 - storage_ratio)^evap_gamma
    @hydroflux flux_name soil_evaporation ~ clamp(potential_evaporation * pet_corr * evap_factor, 0.0, max(0.0, waterstorage))
end

function SOILEVAP_HBV(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    snow_depth::Number=first(@variables snow_depth),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    storage_tension::Number=first(@parameters storage_tension [description = "HBV soil moisture threshold", bounds = (0, 5000), unit = "mm"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_hbv,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(storage_tension, 1.0e-12), 0.0, 1.0)
    @hydroflux flux_name soil_evaporation ~ clamp(
        potential_evaporation * pet_corr * ifelse(snow_depth > 0.0, 0.0, storage_ratio),
        0.0,
        max(0.0, waterstorage),
    )
end

function SOILEVAP_HBV_ORESUND(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    field_capacity::Number=first(@parameters field_capacity [description = "HBV field capacity", bounds = (0, 5000), unit = "mm"]),
    lp::Number=first(@parameters lp [description = "HBV evaporation threshold fraction", bounds = (0, 1), unit = "-"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_hbv_oresund,
)
    threshold_storage = max(lp * field_capacity, 1.0e-12)
    storage_ratio = clamp(max(0.0, waterstorage) / threshold_storage, 0.0, 1.0)
    @hydroflux flux_name soil_evaporation ~ clamp(potential_evaporation * pet_corr * storage_ratio, 0.0, max(0.0, waterstorage))
end

function SOILEVAP_PRMS(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    recharge_storage::Number=first(@variables recharge_storage),
    soil_storage::Number=first(@variables soil_storage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    recharge_max::Number=first(@parameters recharge_max [description = "Maximum recharge-zone storage", bounds = (0, 5000), unit = "mm"]),
    soil_max::Number=first(@parameters soil_max [description = "Maximum lower soil storage", bounds = (0, 5000), unit = "mm"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_prms,
)
    remaining_pet = max(0.0, potential_evaporation * pet_corr)
    upper_ratio = clamp(max(0.0, recharge_storage) / max(recharge_max, 1.0e-12), 0.0, 1.0)
    upper_evap = min(max(0.0, recharge_storage), remaining_pet * upper_ratio)
    lower_pet = max(0.0, remaining_pet - upper_evap)
    lower_ratio = clamp(max(0.0, soil_storage) / max(soil_max, 1.0e-12), 0.0, 1.0)
    lower_evap = ifelse(remaining_pet > max(0.0, recharge_storage), min(max(0.0, soil_storage), lower_pet * lower_ratio), 0.0)
    @hydroflux flux_name soil_evaporation ~ clamp(upper_evap + lower_evap, 0.0, max(0.0, recharge_storage) + max(0.0, soil_storage))
end

function SOILEVAP_SACSMA(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    upper_tension_evaporation::Number=first(@variables upper_tension_evaporation),
    upper_free_evaporation::Number=first(@variables upper_free_evaporation),
    lower_tension_evaporation::Number=first(@variables lower_tension_evaporation),
    uztwc::Number=first(@variables uztwc),
    uzfwc::Number=first(@variables uzfwc),
    lztwc::Number=first(@variables lztwc),
    potential_evaporation::Number=first(@variables potential_evaporation),
    uztwm::Number=first(@parameters uztwm [description = "Upper tension water capacity", bounds = (0, 5000), unit = "mm"]),
    uzfwm::Number=first(@parameters uzfwm [description = "Upper free water capacity", bounds = (0, 5000), unit = "mm"]),
    lztwm::Number=first(@parameters lztwm [description = "Lower tension water capacity", bounds = (0, 5000), unit = "mm"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    upper_tension_flux_name::Symbol=:soilevap_sacsma_uztw,
    upper_free_flux_name::Symbol=:soilevap_sacsma_uzfw,
    lower_tension_flux_name::Symbol=:soilevap_sacsma_lztw,
    total_flux_name::Symbol=:soilevap_sacsma_total,
)
    pet = max(0.0, potential_evaporation * pet_corr)
    euztw = min(max(0.0, uztwc), pet * clamp(max(0.0, uztwc) / max(uztwm, 1.0e-12), 0.0, 1.0))
    euzfw = min(max(0.0, uzfwc), max(0.0, pet - euztw) * clamp(max(0.0, uzfwc) / max(uzfwm, 1.0e-12), 0.0, 1.0))
    elztw = min(max(0.0, lztwc), max(0.0, pet - euztw - euzfw) * clamp(max(0.0, lztwc) / max(lztwm, 1.0e-12), 0.0, 1.0))
    upper_tension_flux = @hydroflux upper_tension_flux_name upper_tension_evaporation ~ euztw
    upper_free_flux = @hydroflux upper_free_flux_name upper_free_evaporation ~ euzfw
    lower_tension_flux = @hydroflux lower_tension_flux_name lower_tension_evaporation ~ elztw
    total_flux = @hydroflux total_flux_name soil_evaporation ~ euztw + euzfw + elztw
    return (upper_tension_flux, upper_free_flux, lower_tension_flux, total_flux)
end

function SOILEVAP_GR4J(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Production store capacity", bounds = (0, 5000), unit = "mm"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_gr4j,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    tanh_term = tanh(max(0.0, potential_evaporation * pet_corr) / max(max_waterstorage, 1.0e-12))
    @hydroflux flux_name soil_evaporation ~ clamp(
        max(0.0, waterstorage) * (2 - storage_ratio) * tanh_term / (1 + (1 - storage_ratio) * tanh_term),
        0.0,
        max(0.0, waterstorage),
    )
end

function SOILEVAP_UBC(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum soil storage", bounds = (0, 5000), unit = "mm"]),
    gamma_e::Number=first(@parameters gamma_e [description = "UBC evaporation parameter", bounds = (1.0e-6, 5000), unit = "mm"]),
    gamma_a::Number=first(@parameters gamma_a [description = "UBC fast-area parameter", bounds = (1.0e-6, 5000), unit = "mm"]),
    imp_fraction::Number=first(@parameters imp_fraction [description = "Impervious fraction", bounds = (0, 1), unit = "-"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_ubc,
)
    beta_fast = imp_fraction * 10^(-(max_waterstorage - waterstorage) / max(gamma_a, 1.0e-12))
    evap_factor = (1 - beta_fast) * 10^(-(max_waterstorage - waterstorage) / max(gamma_e, 1.0e-12))
    @hydroflux flux_name soil_evaporation ~ clamp(max(0.0, potential_evaporation * pet_corr) * max(evap_factor, 0.0), 0.0, max(0.0, waterstorage))
end

function SOILEVAP_PDM(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum soil water storage", bounds = (0, 5000), unit = "mm"]),
    pdm_b::Number=first(@parameters pdm_b [description = "PDM shape parameter", bounds = (0, 20), unit = "-"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_pdm,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    evap_factor = 1 - (1 - storage_ratio)^pdm_b
    @hydroflux flux_name soil_evaporation ~ clamp(max(0.0, potential_evaporation * pet_corr) * evap_factor, 0.0, max(0.0, waterstorage))
end

function SOILEVAP_HYPR(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum soil water storage", bounds = (0, 5000), unit = "mm"]),
    soil_evap_n::Number=first(@parameters soil_evap_n [description = "HYPR soil evaporation exponent", bounds = (0.01, 20), unit = "-"]),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    flux_name::Symbol=:soilevap_hypr,
)
    storage_ratio = clamp(max(0.0, waterstorage) / max(max_waterstorage, 1.0e-12), 0.0, 1.0)
    evap_factor = 1 + storage_ratio - (1 + storage_ratio^soil_evap_n)^(1 / soil_evap_n)
    @hydroflux flux_name soil_evaporation ~ clamp(potential_evaporation * pet_corr * evap_factor, 0.0, max(0.0, waterstorage))
end

function SOILEVAP_CHU(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    potential_evaporation::Number=first(@variables potential_evaporation),
    chu::Number=first(@variables chu),
    chu_maturity::Number=first(@parameters chu_maturity [description = "Crop heat unit at maturity", bounds = (1.0e-6, 1.0e9), unit = "-"]),
    flux_name::Symbol=:soilevap_chu,
)
    @hydroflux flux_name soil_evaporation ~ max(0.0, potential_evaporation) * clamp(chu / max(chu_maturity, 1.0e-12), 0.0, 1.0)
end

function SOILEVAP_AWBM(;
    evap1::Number=first(@variables evap1),
    evap2::Number=first(@variables evap2),
    evap3::Number=first(@variables evap3),
    potential_evaporation::Number=first(@variables potential_evaporation),
    area1::Number=first(@parameters area1 [description = "Area fraction 1", bounds = (0, 1), unit = "-"]),
    area2::Number=first(@parameters area2 [description = "Area fraction 2", bounds = (0, 1), unit = "-"]),
    area3::Number=first(@parameters area3 [description = "Area fraction 3", bounds = (0, 1), unit = "-"]),
    flux1_name::Symbol=:soilevap_awbm_1,
    flux2_name::Symbol=:soilevap_awbm_2,
    flux3_name::Symbol=:soilevap_awbm_3,
)
    f1 = @hydroflux flux1_name evap1 ~ area1 * max(0.0, potential_evaporation)
    f2 = @hydroflux flux2_name evap2 ~ area2 * max(0.0, potential_evaporation)
    f3 = @hydroflux flux3_name evap3 ~ area3 * max(0.0, potential_evaporation)
    return (f1, f2, f3)
end

function SOILEVAP_HYMOD2(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    potential_evaporation::Number=first(@variables potential_evaporation),
    cstar::Number=first(@variables cstar),
    cmax::Number=first(@parameters cmax [description = "Maximum c-star storage", bounds = (1.0e-6, 5000), unit = "mm"]),
    g::Number=first(@parameters g [description = "HYMOD2 minimum evaporation factor", bounds = (0, 1), unit = "-"]),
    kmax::Number=first(@parameters kmax [description = "Maximum evaporation factor", bounds = (0, 10), unit = "-"]),
    shape_c::Number=first(@parameters shape_c [description = "HYMOD2 shape parameter", bounds = (0.01, 20), unit = "-"]),
    flux_name::Symbol=:soilevap_hymod2,
)
    k = kmax * (g + (1 - g) * (clamp(cstar / max(cmax, 1.0e-12), 0.0, 1.0)^shape_c))
    @hydroflux flux_name soil_evaporation ~ max(0.0, k * potential_evaporation)
end

end



