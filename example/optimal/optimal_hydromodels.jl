using CSV, DataFrames, ComponentArrays
using HydroModels
using DataInterpolations
using OptimizationBBO
using Optimization
using ModelingToolkit
using Distributions
using HydroModelTools
using ProgressMeter
using Dates
using Plots
include("../../src/HydroModelLibrary.jl")

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:10000)
lday_vec = df[ts, "dayl(day)"]
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
pet_vec = @. 29.8 * lday_vec * 24 * 0.611 * exp((17.3 * temp_vec) / (temp_vec + 237.3)) / (temp_vec + 273.2)
flow_vec = df[ts, "flow(mm)"]
warm_up = 365
max_iter = 1000

model = HydroModelLibrary.alpine1.model
param_bounds = getbounds.(get_params(model))
random_param_values = map(param_bounds) do param_bound
    rand(Uniform(param_bound[1], param_bound[2]))
end
init_params = ComponentVector(params=NamedTuple{Tuple(get_param_names(model))}(random_param_values))
ps_axes = getaxes(init_params)
init_states = NamedTuple{Tuple(get_state_names(model))}(zeros(length(get_state_names(model)))) |> ComponentVector

input = (P=prcp_vec, Ep=pet_vec, T=temp_vec)
input_matrix = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(model)]))')
config = (solver=HydroModelTools.ODESolver(), interp=LinearInterpolation)
run_kwargs = (config=config, initstates=init_states)
y_mean = mean(flow_vec)
loss_func(y, y_hat) = sum((y .- y_hat) .^ 2) ./ sum((y .- y_mean) .^ 2)

function obj_func(p, _)
    return loss_func(
        flow_vec[warm_up:end],
        model(input_matrix, ComponentVector(p, ps_axes); run_kwargs...)[end, warm_up:end]
    )
end

max_iter = 1000
progress = Progress(max_iter, desc="Optimization")
recorder = []
callback_func!(state, l) = begin
    push!(recorder, (iter=state.iter, loss=l, time=now(), params=state.u))
    next!(progress)
    false
end

lb_list = first.(param_bounds) .|> eltype(input_matrix)
ub_list = last.(param_bounds) .|> eltype(input_matrix)

optf = Optimization.OptimizationFunction(obj_func)
optprob = Optimization.OptimizationProblem(optf, Vector(init_params), lb=lb_list, ub=ub_list)
sol = Optimization.solve(
    optprob,
    BBO_adaptive_de_rand_1_bin_radiuslimited(),
    maxiters=max_iter,
    callback=callback_func!
)

flow_pred = model(input_matrix, ComponentVector(sol.u, ps_axes); run_kwargs...)[end, :]
plot(flow_vec, label="Observed")
plot!(flow_pred, label="Predicted")