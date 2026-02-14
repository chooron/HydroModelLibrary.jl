#! 需要重复率定模型参数看看是否会陷入局部最小

using Pkg

Pkg.activate(".")
using CSV, DataFrames, Dates, ComponentArrays, NPZ
using HydroModels
using DataInterpolations
using Optimization
using OptimizationBBO
using Distributions
using ProgressMeter
using Plots
using JLD2
using Symbolics
using Statistics
using DotEnv

DotEnv.load!()
include("../src/HydroModelLibrary.jl")

# load data
data = npzread(ENV["CAMESL_DATASET_PATH"])
n_catchments = size(data["forcing"], 2)

# build full date index and split ranges
full_dates = collect(Date(1980, 10, 1):Day(1):Date(2014, 9, 30))

train_start = Date(1989, 1, 1)
train_end = Date(1998, 12, 31)
test_start = Date(1999, 1, 1)
test_end = Date(2009, 12, 31)

# align dates to available data length (target may be shorter than full range)
nt = size(data["forcing"], 1)
dates = full_dates[1:nt]
@assert nt == size(data["target"], 1) "Forcing/target length mismatch"

train_idxs = findall(d -> d >= train_start && d <= train_end, dates)
test_idxs = findall(d -> d >= test_start && d <= test_end, dates)
warm_up = 365

@assert length(train_idxs) > warm_up "Training split too short after warm-up"
@assert length(test_idxs) > warm_up "Test split too short after warm-up"

max_iter = 10000
model_nm = "collie1"
const TARGET_BASIN_IDX = 7
const REPEAT_COUNT = 10
@assert 1 <= TARGET_BASIN_IDX <= n_catchments "TARGET_BASIN_IDX out of range"

# build model and reusable metadata once
model_module = HydroModelLibrary.load_model(Symbol(model_nm))
model = model_module.model
model_parameters = model_module.model_parameters
model_params_names = tosymbol.(model_parameters)
model_input_names = HydroModels.get_input_names(model)
model_state_names = HydroModels.get_state_names(model)
model_output_names = HydroModels.get_output_names(model)
model_state_output_names = vcat(model_state_names, model_output_names) |> Tuple
model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModels.getbounds.(model_parameters))

# bounds and config shared by all basins
lb_list = first.(model_params_bounds |> collect)
ub_list = last.(model_params_bounds |> collect)
config = (solver=HydroModels.ODESolver, interp=LinearInterpolation)
result_dir = joinpath(@__DIR__, "calibration_results", "$(model_nm)_repeat", "basin_$(TARGET_BASIN_IDX)")
mkpath(result_dir)
metrics_df = DataFrame(run=Int[], train_kge=Float64[], test_kge=Float64[], train_loss=Float64[], test_loss=Float64[])

kge_score(y, y_hat) = begin
    r = cor(y, y_hat)
    alpha = std(y_hat) / std(y)
    beta = mean(y_hat) / mean(y)
    1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
end

function kge_loss_and_score(y, sim, idxs)
    @assert length(idxs) > warm_up "Not enough points for warm-up in this split"
    eval_idxs = idxs[(warm_up+1):end]
    y_sub = y[eval_idxs]
    sim_sub = sim[end, eval_idxs]
    kge = kge_score(y_sub, sim_sub)
    return 1 - kge, kge
end

split_loss_and_kge(y, sim, idxs) = kge_loss_and_score(y, sim, idxs)

tidx = dayofyear.(dates)
basin_forcing = data["forcing"][:, TARGET_BASIN_IDX, :]
basin_target = data["target"][:, TARGET_BASIN_IDX]
input = (P=basin_forcing[:, 1], T=basin_forcing[:, 2], Ep=basin_forcing[:, 3], Tidx=tidx)
flow_vec = basin_target
input_matrix = stack(input[model_input_names], dims=1)

for run_idx in 1:REPEAT_COUNT
    result_file = joinpath(result_dir, "basin_$(TARGET_BASIN_IDX)_run_$(run_idx).jld2")

    random_param_values = map(zip(model_params_names, model_params_bounds)) do (param_name, param_bound)
        param_name => rand(Uniform(param_bound[1], param_bound[2]))
    end |> NamedTuple
    init_params = ComponentVector(params=NamedTuple{Tuple(model_params_names)}(random_param_values))
    ps_axes = getaxes(init_params)

    obj_func(p, _) = begin
        sim = model(input_matrix, ComponentVector(p, ps_axes); config=config)
        train_loss, _ = split_loss_and_kge(flow_vec, sim, train_idxs)
        train_loss
    end

    progress = Progress(max_iter, desc="Optimization run $(run_idx)/$(REPEAT_COUNT) for basin $(TARGET_BASIN_IDX)")
    recorder = Vector{Any}()
    last_time = time()
    callback_func!(state, l) = begin
        current_time = time()
        delta_time = current_time - last_time
        params_vector = collect(state.u)
        push!(recorder, (iter=state.iter, loss=l, time=now(), delta_time=delta_time, params=params_vector))
        last_time = current_time
        next!(progress)
        false
    end

    optf = Optimization.OptimizationFunction(obj_func, Optimization.AutoForwardDiff())
    optprob = Optimization.OptimizationProblem(optf, Vector(init_params), lb=lb_list, ub=ub_list)
    start_time = time()
    sol = Optimization.solve(
        optprob,
        BBO_adaptive_de_rand_1_bin_radiuslimited(),
        maxiters=max_iter,
        callback=callback_func!
    )
    elapsed_seconds = round(Int, time() - start_time)

    params = ComponentVector(sol.u, ps_axes)
    output = model(input_matrix, params; config=config)

    train_loss, train_kge = split_loss_and_kge(flow_vec, output, train_idxs)
    test_loss, test_kge = split_loss_and_kge(flow_vec, output, test_idxs)

    recorder_df = DataFrame(recorder)
    @assert !isempty(recorder_df) "No optimization steps recorded"
    loss_history = recorder_df.loss
    cumulative_time = cumsum(recorder_df.delta_time)

    param_matrix = reduce(hcat, recorder_df.params)
    param_df = DataFrame(permutedims(param_matrix), model_params_names)
    trajectory_df = hcat(select(recorder_df, Not(:params)), param_df)

    @save result_file params loss_history train_kge test_kge train_loss test_loss elapsed_seconds cumulative_time recorder_df trajectory_df

    push!(metrics_df, (run=run_idx, train_kge=train_kge, test_kge=test_kge, train_loss=train_loss, test_loss=test_loss))
    @info "run metrics" basin=TARGET_BASIN_IDX run=run_idx train_kge=train_kge test_kge=test_kge train_loss=train_loss test_loss=test_loss

    GC.gc()
end

@info "repeat metrics" metrics_df
batch_metrics_file = joinpath(result_dir, "repeat_metrics.jld2")
@save batch_metrics_file metrics_df