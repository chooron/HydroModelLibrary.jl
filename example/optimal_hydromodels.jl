using CSV, DataFrames, Dates, ComponentArrays
using HydroModels, HydroModelTools
using DataInterpolations
using Optimization
using OptimizationBBO
using Distributions
using ProgressMeter
using Plots
using JLD2
using Symbolics
include("../src/HydroModelLibrary.jl")

# load data
df = CSV.read("data/marrmot/3604000.csv", DataFrame)
tidx = rem.(0:length(df[!, "prec"])-1, 365) .+ 1
input = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"], Tidx=tidx)
flow_vec = df[!, "flow"]
warm_up = 365
max_iter = 10000
result_dict = Dict{Symbol, Any}()

for model_nm in HydroModelLibrary.AVAILABLE_MODELS
    if isdir("cache/calibrate/$model_nm")
        @info "cache/calibrate/$model_nm already exists, skip calibration"
        continue
    else
        mkdir("cache/calibrate/$model_nm")
    end

    model_module = HydroModelLibrary.load_model(Symbol(model_nm))
    model = model_module.model
    model_parameters = model_module.model_parameters
    model_params_names = tosymbol.(model_parameters)
    model_input_names = HydroModels.get_input_names(model)
    model_state_names = HydroModels.get_state_names(model)
    model_output_names = HydroModels.get_output_names(model)

    model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModels.getbounds.(model_parameters))
    random_param_values = map(zip(model_params_names, model_params_bounds)) do (param_name, param_bound)
        param_name => rand(Uniform(param_bound[1], param_bound[2]))
    end |> NamedTuple
    init_params = ComponentVector(params=NamedTuple{Tuple(model_params_names)}(random_param_values))
    ps_axes = getaxes(init_params)

    input_matrix = stack(input[model_input_names], dims=1)
    config = (solver=HydroModelTools.ODESolver(), interp=LinearInterpolation)

    y_mean = mean(flow_vec)
    mse_func(y, y_hat) = sum((y .- y_hat) .^ 2) ./ length(y)
    r2_func(y, y_hat) = sum((y .- y_hat) .^ 2) ./ sum((y .- y_mean) .^ 2)
    kge_func(y, y_hat) = sqrt((r2_func(y, y_hat))^2 + (std(y_hat) / std(y) - 1)^2 + (mean(y_hat) / mean(y) - 1)^2)
    function obj_func(p, _)
        return r2_func(
            flow_vec[warm_up:end],
            model(input_matrix, ComponentVector(p, ps_axes); config=config)[end, warm_up:end]
        )
    end

    progress = Progress(max_iter, desc="Optimization")
    recorder = []
    callback_func!(state, l) = begin
        push!(recorder, (iter=state.iter, loss=l, time=now(), params=state.u))
        next!(progress)
        false
    end

    lb_list = first.(model_params_bounds |> collect) .|> eltype(input_matrix)
    ub_list = last.(model_params_bounds |> collect) .|> eltype(input_matrix)

    optf = Optimization.OptimizationFunction(obj_func, Optimization.AutoForwardDiff())
    optprob = Optimization.OptimizationProblem(optf, Vector(init_params), lb=lb_list, ub=ub_list)
    sol = Optimization.solve(
        optprob,
        BBO_adaptive_de_rand_1_bin_radiuslimited(),
        maxiters=max_iter,
        callback=callback_func!
    )
    recorder_df = DataFrame(recorder)
    CSV.write("cache/calibrate/$model_nm/recorder.csv", recorder_df)
    params = ComponentVector(sol.u, ps_axes)
    JLD2.save("cache/calibrate/$model_nm/sol.jld2", "opt_params", params)
    output = model(input_matrix, params; config=config)
    model_output_names = vcat(model_state_names, model_output_names) |> Tuple
    output_df = DataFrame(NamedTuple{model_output_names}(eachslice(output, dims=1)))
    CSV.write("cache/calibrate/$model_nm/output.csv", output_df)
    loss = r2_func(flow_vec, output_df[!, :Qt])
    result_dict[model_nm] = loss
end