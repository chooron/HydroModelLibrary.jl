module HydroModelLibrary

# using SpecialFunctions

using HydroModels
using HydroModelCore

using ComponentArrays

using NNlib
using Symbolics: tosymbol

## -------------------------- smooth function -------------------------- ##
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
smoothlogistic_func(S, Smax, r=0.01, e=5.0) = 1 / (1 + exp((S - r * e * Smax) / (r * Smax)))

## -------------------------- hydrological model -------------------------- ##

function load_model(model_name::Symbol; reload=false)
    if !isdefined(HydroModelLibrary, model_name) || reload
        model_path = joinpath(@__DIR__, "models", "$(model_name).jl")
        include(model_path)
    end
    return getfield(HydroModelLibrary, model_name)
end

function get_params_bounds(model_name::Symbol)
    if !isdefined(HydroModelLibrary, model_name)
        load_model(model_name)
    end
    model_parameters = getfield(HydroModelLibrary, model_name).model_parameters
    model_params_names = tosymbol.(model_parameters)
    model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModelCore.getbounds.(model_parameters))
    return model_params_bounds
end

function get_random_params(model_name::Symbol)
    if !isdefined(HydroModelLibrary, model_name)
        load_model(model_name)
    end
    model_parameters = getfield(HydroModelLibrary, model_name).model_parameters
    model_params_names = tosymbol.(model_parameters)
    model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModelCore.getbounds.(model_parameters))
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
    :echo,
    :exphydro,
    :flexb,
    :gr4j,
    :gsfb,
    :gsmsocont,
    :hbv_edu,
    :hbv,
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

function load_router(router_name::Symbol)
    if router_name in AVAILABLE_ROUTERS
        if !isdefined(HydroModelLibrary, router_name)
            router_path = joinpath(@__DIR__, "routers", "$(router_name).jl")
            include(router_path)
        end
        return getfield(HydroModelLibrary, router_name)
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
