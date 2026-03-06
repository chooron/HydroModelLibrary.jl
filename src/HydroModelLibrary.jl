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

function get_params_bounds(model_name::Symbol)
    model_module = load_model(model_name)
    model_parameters = model_module.model_parameters
    model_params_names = tosymbol.(model_parameters)
    model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModels.getbounds.(model_parameters))
    return model_params_bounds
end

function get_random_params(model_name::Symbol)
    model_module = load_model(model_name)
    model_parameters = model_module.model_parameters
    model_params_names = tosymbol.(model_parameters)
    model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModels.getbounds.(model_parameters))
    random_param_values = map(zip(model_params_names, model_params_bounds)) do (param_name, param_bound)
        param_name => rand(param_bound[1]:param_bound[2])
    end |> NamedTuple
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
    :gamma_uh
]
export AVAILABLE_MODELS, AVAILABLE_ROUTERS

include("routes/gamma_uh.jl")
export GammaHydrograph

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
            router_path = joinpath(@__DIR__, "routers", "$(router_name).jl")
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

