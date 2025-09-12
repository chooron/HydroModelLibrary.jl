module HydroModelLibrary

# using SpecialFunctions

using HydroModels
using HydroModelCore

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
smoothlogistic_func(S, Smax, r=0.01, e=5.0) = 1 / (1 + exp((S - r * e * Smax) / (r * Smax)))

# 定义一个函数来按需加载模型
function load_model(model_name::Symbol; reload=false)
    # 检查模块是否已经加载，如果没有则加载
    if !isdefined(HydroModelLibrary, model_name) || reload
        model_path = joinpath(@__DIR__, "models", "$(model_name).jl")
        include(model_path)
    end
    # 返回模块
    return getfield(HydroModelLibrary, model_name)
end

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

function load_model(model_name::Symbol)
    if model_name in AVAILABLE_MODELS
        if !isdefined(HydroModelLibrary, model_name)
            model_path = joinpath(@__DIR__, "models", "$(model_name).jl")
            include(model_path)
        end
        return getfield(HydroModelLibrary, model_name)
    else
        throw(ArgumentError("Model $model_name is not available"))
    end
end

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

end
