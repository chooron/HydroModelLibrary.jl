using CSV, DataFrames, Dates, ComponentArrays, NPZ
using HydroModels, HydroModelTools
using DataInterpolations
using Optimization
using OptimizationBBO
using Distributions
using ProgressMeter
using Plots
using JLD2
using Symbolics
using DotEnv

DotEnv.load!()
include("../src/HydroModelLibrary.jl")

# load data
data = npzread(ENV["CAMESL_DATASET_PATH"])
gage_ids = data["gage_ids"]
select_idx = findfirst(x -> x == 1013500, gage_ids)
test_forcing = data["forcings"][select_idx, :, :]
area = data["attributes"][select_idx, 12]
test_target = data["target"][select_idx, :]
tidx = rem.(0:size(test_forcing)[1]-1, 365) .+ 1
input = (P=test_forcing[:, 1], T=test_forcing[:, 2], Ep=test_forcing[:, 3], Tidx=tidx)
flow_vec = (10^3) * test_target * 0.0283168 * 3600 * 24 / (area * (10^6))
camels_input = stack(input, dims=1)

warm_up = 365
max_iter = 1000
model_nm = "hbv_edu"

model_module = HydroModelLibrary.load_model(Symbol(model_nm))
model = model_module.model
model_parameters = model_module.model_parameters
model_params_names = tosymbol.(model_parameters)
model_input_names = HydroModels.get_input_names(model)
model_state_names = HydroModels.get_state_names(model)
model_output_names = HydroModels.get_output_names(model)
model_state_output_names = vcat(model_state_names, model_output_names) |> Tuple

model_params_bounds = HydroModelLibrary.get_params_bounds(model_nm)

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
# sol = Optimization.solve(
#     optprob,
#     BBO_adaptive_de_rand_1_bin_radiuslimited(),
#     maxiters=max_iter,
#     callback=callback_func!
# )
# recorder_df = DataFrame(recorder)
# params = ComponentVector(sol.u, ps_axes)
# output = model(input_matrix, params; config=config)
# output_df = DataFrame(NamedTuple{model_state_output_names}(eachslice(output, dims=1)))
# @info "r2" 1 - r2_func(flow_vec, output[end, :])
# plot(output_df[2000:3000, :Qt], label="predicted")
# plot!(flow_vec[2000:3000], label="observed")