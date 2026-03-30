module HydroModelLibrary

# using SpecialFunctions

using HydroModels
using ComponentArrays

using NNlib
using HydroModels: tosymbol

## -------------------------- smooth function -------------------------- ##
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
smoothlogistic_func(S, Smax, r=0.01, e=5.0) = 1 / (1 + exp((S - r * e * Smax) / (r * Smax)))

## -------------------------- hydrological model -------------------------- ##

# Access dynamically included model/router module bindings in the latest world.
_get_library_binding(name::Symbol) = Base.invokelatest(getfield, HydroModelLibrary, name)

function load_model(model_name::Symbol; reload=false)
    if !isdefined(HydroModelLibrary, model_name) || reload
        model_path = joinpath(@__DIR__, "models", "$(model_name).jl")
        include(model_path)
    end
    return _get_library_binding(model_name)
end

function _get_model_parameters(model_name::Symbol)
    model_module = load_model(model_name)

    model_parameters = try
        Base.invokelatest(getproperty, model_module, :model_parameters)
    catch err
        if err isa UndefVarError
            throw(ArgumentError(
                "Model '$model_name' does not define model_parameters in its module wrapper.",
            ))
        end
        rethrow(err)
    end

    model_parameters isa AbstractVector || throw(ArgumentError(
        "Model '$model_name' has invalid model_parameters type: $(typeof(model_parameters)).",
    ))
    isempty(model_parameters) &&
        throw(ArgumentError("Model '$model_name' exposes an empty model_parameters list."))

    return model_parameters
end

function _get_param_names_and_bounds(model_name::Symbol, model_parameters)
    param_names = Symbol.(tosymbol.(model_parameters))
    bounds = Vector{Tuple{Float64,Float64}}(undef, length(model_parameters))

    for (idx, param) in enumerate(model_parameters)
        raw_bounds = try
            Base.invokelatest(HydroModels.getbounds, param)
        catch
            throw(ArgumentError(
                "Missing bounds for parameter '$param' in model '$model_name'.",
            ))
        end

        length(raw_bounds) == 2 || throw(ArgumentError(
            "Invalid bounds for parameter '$param' in model '$model_name': expected 2 values.",
        ))

        lo = Float64(raw_bounds[1])
        hi = Float64(raw_bounds[2])
        lo <= hi || throw(ArgumentError(
            "Invalid bounds for parameter '$param' in model '$model_name': lower bound $lo is greater than upper bound $hi.",
        ))

        bounds[idx] = (lo, hi)
    end

    return param_names, bounds
end

function get_params_bounds(model_name::Symbol)
    model_parameters = _get_model_parameters(model_name)
    param_names, bounds = _get_param_names_and_bounds(model_name, model_parameters)
    return NamedTuple{Tuple(param_names)}(Tuple(bounds))
end

function get_random_params(model_name::Symbol)
    model_parameters = _get_model_parameters(model_name)
    param_names, bounds = _get_param_names_and_bounds(model_name, model_parameters)

    sampled_values = ntuple(length(param_names)) do i
        lo, hi = bounds[i]
        lo == hi ? lo : lo + (hi - lo) * rand()
    end

    random_param_values = NamedTuple{Tuple(param_names)}(sampled_values)
    init_params = ComponentVector(params=random_param_values)
    return init_params
end

load_model(model_name::String; reload=false) = load_model(Symbol(model_name); reload=reload)
get_random_params(model_name::String) = get_random_params(Symbol(model_name))
get_params_bounds(model_name::String) = get_params_bounds(Symbol(model_name))

AVAILABLE_MODELS = [
    :alpine1,
    :alpine2,
    :australia,
    :collie1,
    :collie2,
    :collie3,
    :cemaneigegr4j,
    :echo,
    :exphydro,
    :flexb,
    :gr4j,
    :gsfb,
    :gsmsocont,
    :hbv_edu,
    :hbv,
    :hbv96,
    :hillslope,
    :hmets,
    :hymod,
    :ihacres,
    :ihm19,
    :lascam,
    :mcrm,
    :modhydrolog,
    :mopex1,
    :mopex2,
    :mopex3,
    :mopex4,
    :mopex5,
    :nam,
    :newzealand1,
    :newzealand2,
    :penman,
    :plateau,
    :prms,
    :sacramento,
    # :smar, nash cascade routing is not implemented
    # :smart, skip this model
    :susannah1,
    :susannah2,
    :tank,
    :tcm,
    :unitedstates,
    :wetland,
    :xaj,
    # :xinanjiang # bad model
]

AVAILABLE_ROUTERS = [
    :rapid,
    :gamma_uh,
    :unithydro,
    :channel_route,
]
export AVAILABLE_MODELS, AVAILABLE_ROUTERS

include("routes/unithydro.jl")
include("routes/channel_route.jl")
include("routes/gamma_uh.jl")
using .gamma_uh: GammaHydrograph
using .unithydro: ROUTE_DUMP,
                  ROUTE_GAMMA_CONVOLUTION,
                  ROUTE_TRI_CONVOLUTION,
                  ROUTE_RESERVOIR_SERIES,
                  ROUTE_DIFFUSIVE_WAVE,
                  ROUTE_NONE,
                  ROUTE_PLUG_FLOW,
                  ROUTE_LINEAR_STORAGE,
                  ROUTE_STORAGE_COEFF,
                  ROUTE_MUSKINGUM,
                  ROUTE_NASH_CASCADE,
                  ROUTE_HYDROLOGIC,
                  TOC_MCDERMOTT_PILGRIM,
                  TOC_BRANSBY_WILLIAMS,
                  TOC_WILLIAMS_1922,
                  build_unit_hydrograph,
                  dump_unit_hydrograph,
                  gamma_unit_hydrograph,
                  triangular_unit_hydrograph,
                  nash_unit_hydrograph,
                  toc_mcdermott_pilgrim,
                  toc_bransby_williams,
                  toc_williams_1922,
                  time_of_concentration
using .channel_route: ChannelRoute,
                      build_channel_route,
                      no_channel_route,
                      plug_flow_route,
                      diffusive_wave_route,
                      linear_storage_route,
                      storage_coeff_route,
                      muskingum_route,
                      nash_cascade_route,
                      hydrologic_route
export @channelroute,
       GammaHydrograph,
       ROUTE_DUMP,
       ROUTE_GAMMA_CONVOLUTION,
       ROUTE_TRI_CONVOLUTION,
       ROUTE_RESERVOIR_SERIES,
       ROUTE_DIFFUSIVE_WAVE,
       ROUTE_NONE,
       ROUTE_PLUG_FLOW,
       ROUTE_LINEAR_STORAGE,
       ROUTE_STORAGE_COEFF,
       ROUTE_MUSKINGUM,
       ROUTE_NASH_CASCADE,
       ROUTE_HYDROLOGIC,
       TOC_MCDERMOTT_PILGRIM,
       TOC_BRANSBY_WILLIAMS,
       TOC_WILLIAMS_1922,
       build_unit_hydrograph,
       dump_unit_hydrograph,
       gamma_unit_hydrograph,
       triangular_unit_hydrograph,
       nash_unit_hydrograph,
       toc_mcdermott_pilgrim,
       toc_bransby_williams,
       toc_williams_1922,
       time_of_concentration,
       ChannelRoute,
       build_channel_route,
       no_channel_route,
       plug_flow_route,
       diffusive_wave_route,
       linear_storage_route,
       storage_coeff_route,
       muskingum_route,
       nash_cascade_route,
       hydrologic_route

## -------------------------- hydrological fluxes -------------------------- ##
include("fluxes/rainsnow.jl")
include("fluxes/precipinterception.jl")
include("fluxes/abstraction.jl")
include("fluxes/baseflow.jl")
include("fluxes/infiltration.jl")
include("fluxes/capillaryrise.jl")
include("fluxes/soilbalance.jl")
include("fluxes/canopyevap.jl")
include("fluxes/canopydrip.jl")
include("fluxes/openwaterevap.jl")
include("fluxes/overflow.jl")
include("fluxes/depressionstorage.jl")
include("fluxes/seepage.jl")
include("fluxes/lakerelease.jl")
include("fluxes/percolation.jl")
include("fluxes/potentialmelt.jl")
include("fluxes/interflow.jl")
include("fluxes/bottomdrain.jl")
include("fluxes/snowmelt.jl")
include("fluxes/snowsublimation.jl")
include("fluxes/snowfreeze.jl")
include("fluxes/snowalbedo.jl")
include("fluxes/snowbalance.jl")
include("fluxes/glaciermelt.jl")
include("fluxes/glacierrelease.jl")
include("fluxes/glacierinfiltration.jl")
include("fluxes/lakefreeze.jl")
include("fluxes/cropheatunit.jl")
include("fluxes/specialprocesses.jl")
include("fluxes/soilevap.jl")

using .RainSnow: RAINSNOW_DATA, RAINSNOW_DINGMAN, RAINSNOW_HBV, RAINSNOW_UBCWM
using .PrecipInterception: PRECIP_ICEPT_LAI, PRECIP_ICEPT_EXPLAI
using .Abstraction: ABST_FILL, ABST_MAX, ABST_RATIO, ABST_DERIVED, ABST_SCS, ABST_PERCENTAGE, ABST_PDMROF, ABST_UWFS
using .Baseflow: BASE_CONSTANT,
                 BASE_LINEAR,
                 BASE_LINEAR_ANALYTIC,
                 BASE_POWER_LAW,
                 BASE_VIC,
                 BASE_GR4J,
                 BASE_TOPMODEL,
                 BASE_THRESH_POWER,
                 BASE_THRESH_STOR
using .Infiltration: INF_RATIONAL,
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
using .CapillaryRise: CAPRISE_HBV, CRISE_HBV
using .SoilBalance: SOILBAL_BUCKET, SOILBAL_PRMS, SOILBAL_PENMAN
using .CanopyEvap: CANEVAP_ALL, CANEVAP_LINEAR, CANEVAP_PRMS, CANEVP_MAXIMUM, CANEVP_RUTTER
using .CanopyDrip: CANDRIP_EXCESS, CANDRIP_PRMS, CANDRIP_LINEAR, CANDRIP_SLOWDRAIN
using .OpenWaterEvap: OWEVAP_DATA,
                      OWEVAP_HARGREAVES_1985,
                      OWEVAP_HARGREAVES,
                      OWEVAP_PENMAN_SIMPLE,
                      OWEVAP_PENMAN_MONTEITH,
                      OPEN_WATER_EVAP,
                      OPEN_WATER_RIPARIAN,
                      OPEN_WATER_UWFS
using .Overflow: OVERFLOW_ALL, OVERFLOW_THRESHOLD, OVERFLOW_LINEAR, OVERFLOW_NONLINEAR, OVERFLOW_GR4J
using .DepressionStorage: DEPSTOR_BUCKET, DEPSTOR_PRMS, DFLOW_THRESHPOW, DFLOW_LINEAR, DFLOW_WEIR
using .Seepage: SEEPAGE_LINEAR, SEEP_LINEAR, SEEPAGE_THRESHOLD
using .LakeRelease: LAKEREL_LINEAR
using .Percolation: PERC_CONSTANT,
                    PERC_LINEAR,
                    PERC_POWER_LAW,
                    PERC_PRMS,
                    PERC_SACRAMENTO,
                    PERC_GAWSER,
                    PERC_GAWSER_CONSTRAIN,
                    PERC_GR4JEXCH,
                    PERC_GR4JEXCH2
using .PotentialMelt: POTMELT_DATA, POTMELT_DEGREE_DAY, POTMELT_DD_FREEZE, POTMELT_HBV, POTMELT_HBV_ROS, POTMELT_HMETS, POTMELT_RESTRICTED
using .Interflow: INTERFLOW_PRMS
using .BottomDrain: BOTTOMDRAIN_LINEAR, BOTTOMDRAIN_POWER, BOTTOMDRAIN_THRESH
using .SnowMelt: SNOWMELT_DEGREE_DAY, SNOWMELT_FROM_POTENTIAL
using .SnowSublimation: SUBLIM_KUZMIN, SUBLIM_CENTRAL_SIERRA, SUBLIM_BULK_AERO
using .SnowFreeze: FREEZE_DEGREE_DAY, FREEZE_HMETS
using .SnowAlbedo: SNOALB_UBC, SNOALB_CRHM_ESSERY, SNOALB_BAKER
using .SnowBalance: SNOBAL_SIMPLE_MELT, SNOBAL_HBV, SNOBAL_HMETS, SNOBAL_CEMA_NIEGE, SNOBAL_CEMA_NEIGE, SNOBAL_COLD_CONTENT, SNOBAL_TWO_LAYER
using .GlacierMelt: GLACIERMELT_DEGREE_DAY, GLACIERMELT_FROM_POTENTIAL, GLACIERMELT_GSMSOCONT
using .GlacierRelease: GRELEASE_LINEAR_STORAGE, GRELEASE_LINEAR_ANALYTIC, GRELEASE_HBV_EC
using .GlacierInfiltration: GLACIERINF_ALL, GLACIERINF_LINEAR, GLACIERINF_GSMSOCONT
using .LakeFreeze: LFREEZE_BASIC
using .CropHeatUnit: CHU_ONTARIO
using .SpecialProcesses: CONVOL_GR4J_1, CONVOL_GR4J_2, CONVOL_GAMMA, CONVOL_GAMMA2
using .SoilEvap: SOILEVAP_ALL,
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

export RAINSNOW_DATA,
       RAINSNOW_DINGMAN,
       RAINSNOW_HBV,
       RAINSNOW_UBCWM,
       PRECIP_ICEPT_LAI,
       PRECIP_ICEPT_EXPLAI,
       ABST_FILL,
       ABST_MAX,
       ABST_RATIO,
       ABST_DERIVED,
       ABST_SCS,
       ABST_PERCENTAGE,
       ABST_PDMROF,
       ABST_UWFS,
       BASE_CONSTANT,
       BASE_LINEAR,
       BASE_LINEAR_ANALYTIC,
       BASE_POWER_LAW,
       BASE_VIC,
       BASE_GR4J,
       BASE_TOPMODEL,
       BASE_THRESH_POWER,
       BASE_THRESH_STOR,
       INF_RATIONAL,
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
       INF_PDM,
       CAPRISE_HBV,
       CRISE_HBV,
       SOILBAL_BUCKET,
       SOILBAL_PRMS,
       SOILBAL_PENMAN,
       CANEVAP_ALL,
       CANEVAP_LINEAR,
       CANEVAP_PRMS,
       CANEVP_MAXIMUM,
       CANEVP_RUTTER,
       CANDRIP_EXCESS,
       CANDRIP_PRMS,
       CANDRIP_LINEAR,
       CANDRIP_SLOWDRAIN,
       OWEVAP_DATA,
       OWEVAP_HARGREAVES_1985,
       OWEVAP_HARGREAVES,
       OWEVAP_PENMAN_SIMPLE,
       OWEVAP_PENMAN_MONTEITH,
       OPEN_WATER_EVAP,
       OPEN_WATER_RIPARIAN,
       OPEN_WATER_UWFS,
       OVERFLOW_ALL,
       OVERFLOW_THRESHOLD,
       OVERFLOW_LINEAR,
       OVERFLOW_NONLINEAR,
       OVERFLOW_GR4J,
       DEPSTOR_BUCKET,
       DEPSTOR_PRMS,
       DFLOW_THRESHPOW,
       DFLOW_LINEAR,
       DFLOW_WEIR,
       SEEPAGE_LINEAR,
       SEEP_LINEAR,
       SEEPAGE_THRESHOLD,
       LAKEREL_LINEAR,
       PERC_CONSTANT,
       PERC_LINEAR,
       PERC_POWER_LAW,
       PERC_PRMS,
       PERC_SACRAMENTO,
       PERC_GAWSER,
       PERC_GAWSER_CONSTRAIN,
       PERC_GR4JEXCH,
       PERC_GR4JEXCH2,
       POTMELT_DATA,
       POTMELT_DEGREE_DAY,
       POTMELT_DD_FREEZE,
       POTMELT_HBV,
       POTMELT_HBV_ROS,
       POTMELT_HMETS,
       POTMELT_RESTRICTED,
       INTERFLOW_PRMS,
       BOTTOMDRAIN_LINEAR,
       BOTTOMDRAIN_POWER,
       BOTTOMDRAIN_THRESH,
       SNOWMELT_DEGREE_DAY,
       SNOWMELT_FROM_POTENTIAL,
       SUBLIM_KUZMIN,
       SUBLIM_CENTRAL_SIERRA,
       SUBLIM_BULK_AERO,
       FREEZE_DEGREE_DAY,
       FREEZE_HMETS,
       SNOALB_UBC,
       SNOALB_CRHM_ESSERY,
       SNOALB_BAKER,
       SNOBAL_SIMPLE_MELT,
       SNOBAL_HBV,
       SNOBAL_HMETS,
       SNOBAL_CEMA_NIEGE,
       SNOBAL_CEMA_NEIGE,
       SNOBAL_COLD_CONTENT,
       SNOBAL_TWO_LAYER,
       GLACIERMELT_DEGREE_DAY,
       GLACIERMELT_FROM_POTENTIAL,
       GLACIERMELT_GSMSOCONT,
       GRELEASE_LINEAR_STORAGE,
       GRELEASE_LINEAR_ANALYTIC,
       GRELEASE_HBV_EC,
       GLACIERINF_ALL,
       GLACIERINF_LINEAR,
       GLACIERINF_GSMSOCONT,
       LFREEZE_BASIC,
       CHU_ONTARIO,
       CONVOL_GR4J_1,
       CONVOL_GR4J_2,
       CONVOL_GAMMA,
       CONVOL_GAMMA2,
       SOILEVAP_ALL,
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

## -------------------------- raven template models -------------------------- ##
include("raven-model/hbv_light.jl")
include("raven-model/gr4j.jl")
include("raven-model/hymod.jl")
include("raven-model/todo_models.jl")

using .RavenHBVLight: build_raven_hbv_light
using .RavenGR4J: build_raven_gr4j
using .RavenHYMOD: build_raven_hymod
using .RavenTemplateTODOs: RAVEN_TEMPLATE_MODEL_STATUS,
                           build_raven_ubcwm,
                           build_raven_hbv_ec,
                           build_raven_canadian_shield,
                           build_raven_mohyse,
                           build_raven_hmets,
                           build_raven_hypr,
                           build_raven_awbm,
                           build_raven_sac_sma,
                           build_raven_routing_only,
                           build_raven_blended_v1,
                           build_raven_blended_v2

AVAILABLE_RAVEN_MODELS = [
    :ubcwm,
    :hbv_ec,
    :hbv_light,
    :gr4j,
    :canadian_shield,
    :mohyse,
    :hmets,
    :hypr,
    :hymod,
    :awbm,
    :sac_sma,
    :routing_only,
    :blended_v1,
    :blended_v2,
]

export AVAILABLE_RAVEN_MODELS,
       RAVEN_TEMPLATE_MODEL_STATUS,
       build_raven_hbv_light,
       build_raven_gr4j,
       build_raven_hymod,
       build_raven_ubcwm,
       build_raven_hbv_ec,
       build_raven_canadian_shield,
       build_raven_mohyse,
       build_raven_hmets,
       build_raven_hypr,
       build_raven_awbm,
       build_raven_sac_sma,
       build_raven_routing_only,
       build_raven_blended_v1,
       build_raven_blended_v2

## -------------------------- adequacy assessment -------------------------- ##
include("criteria/adequacy.jl")
export ADEQUACY_PE_LIMIT,
       percentage_error, is_adequate,
       monthly_mean_daily_flow, yearly_mean_daily_flow,
       runoff_ratio, bfi_ioh1980, mhfd, slope_fdc,
       evaluate_adequacy_assessment

## -------------------------- accuracy assessment -------------------------- ##
include("criteria/accuracy.jl")
export acc, d, d1, d1_p, dmod, dr, drel, ed, g_mean_diff,
       h1_mhe, h1_mahe, h1_rmshe,
       h10_mhe, h10_mahe, h10_rmshe,
       h2_mhe, h2_mahe, h2_rmshe,
       h3_mhe, h3_mahe, h3_rmshe,
       h4_mhe, h4_mahe, h4_rmshe,
       h5_mhe, h5_mahe, h5_rmshe,
       h6_mhe, h6_mahe, h6_rmshe,
       h7_mhe, h7_mahe, h7_rmshe,
       h8_mhe, h8_mahe, h8_rmshe,
       irmse, kge_2009, kge_2012, lm_index,
       maappe, mae, male, mapd, mape, mase, mdae, mde, mdse, me, mean_var,
       mhe_h1, mhe_h2, mhe_h3, mhe_h4,
       mle, mse, msle, ned,
       nrmse_range, nrmse_mean, nrmse_iqr,
       nse, nse_mod, nse_rel,
       person_r, spearman_r, r_squared, mb_r, r2,
       rmse, rmsle, sa, sc, sga, sid, smape1, smape2, ve, watt_m

function load_router(router_name::Symbol)
    if router_name in AVAILABLE_ROUTERS
        if !isdefined(HydroModelLibrary, router_name)
            router_path = joinpath(@__DIR__, "routes", "$(router_name).jl")
            include(router_path)
        end
        return _get_library_binding(router_name)
    else
        throw(ArgumentError("Router $router_name is not available"))
    end
end
export load_model, load_router

## -------------------------- sample data loading -------------------------- ##
include("sample_data.jl")
export load_sample_data, list_sample_data, get_sample_data_info, 
       load_sample_data_for_model, AVAILABLE_SAMPLE_DATA

end































