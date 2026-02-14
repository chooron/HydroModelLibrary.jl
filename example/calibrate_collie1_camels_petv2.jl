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
const USE_INVERSE_KGE = true
const INVKGE_EPS_FACTOR = 0.01
const INVKGE_MIN_EPS = 1e-3
const INVKGE_STABILITY_EPS = 1e-5

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
result_dir = joinpath(@__DIR__, "calibration_results", "collie1_invkge")
mkpath(result_dir)

metrics_df = DataFrame(basin=Int[], train_kge=Float64[], test_kge=Float64[], train_loss=Float64[], test_loss=Float64[])
all_recorders = Vector{Any}(undef, n_catchments)
all_params = Vector{Any}(undef, n_catchments)

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

function inverse_kge_loss_and_score(
    y,
    sim,
    idxs;
    eps_factor = INVKGE_EPS_FACTOR,
    min_eps = INVKGE_MIN_EPS,
    stability_eps = INVKGE_STABILITY_EPS,
)
    @assert length(idxs) > warm_up "Not enough points for warm-up in this split"
    eval_idxs = idxs[(warm_up+1):end]
    y_sub = y[eval_idxs]
    sim_sub = sim[end, eval_idxs]

    obs_mean = mean(y_sub)
    eps = max(eps_factor * obs_mean, min_eps)

    y_safe = max.(y_sub, 0.0)
    sim_safe = max.(sim_sub, 0.0)

    y_inv = 1 ./(y_safe .+ eps)
    sim_inv = 1 ./(sim_safe .+ eps)

    mean_t = mean(y_inv)
    mean_p = mean(sim_inv)
    std_t = std(y_inv)
    std_p = std(sim_inv)

    dev_t = y_inv .- mean_t
    dev_p = sim_inv .- mean_p
    numerator = sum(dev_p .* dev_t)
    denominator = sqrt(sum(dev_p .^ 2) * sum(dev_t .^ 2)) + stability_eps
    r = numerator / denominator

    beta = mean_p / (mean_t + stability_eps)
    gamma = std_p / (std_t + stability_eps)

    kge = 1 - sqrt((r - 1)^2 + (beta - 1)^2 + (gamma - 1)^2)
    return 1 - kge, kge
end

split_loss_and_kge(y, sim, idxs) = USE_INVERSE_KGE ? inverse_kge_loss_and_score(y, sim, idxs) : kge_loss_and_score(y, sim, idxs)



for select_idx in 1:n_catchments
    result_file = joinpath(result_dir, "basin_$(select_idx).jld2")

    # 如果已存在结果文件，直接加载并记录指标，不再训练
    if isfile(result_file)
        saved = JLD2.load(result_file)
        push!(metrics_df, (
            basin=select_idx,
            train_kge=get(saved, :train_kge, NaN),
            test_kge=get(saved, :test_kge, NaN),
            train_loss=get(saved, :train_loss, NaN),
            test_loss=get(saved, :test_loss, NaN)
        ))
        all_params[select_idx] = get(saved, :params, nothing)
        all_recorders[select_idx] = nothing
        @info "basin already processed; skip" basin = select_idx
    else
        # basin-specific data prep
        test_forcing = data["forcing"][:, select_idx, :]
        test_target = data["target"][:, select_idx]
        tidx = dayofyear.(dates)
        input = (P=test_forcing[:, 1], T=test_forcing[:, 2], Ep=test_forcing[:, 3], Tidx=tidx)
        flow_vec = test_target

        # fresh random starting point per basin
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

        input_matrix = stack(input[model_input_names], dims=1)

        progress = Progress(max_iter, desc="Optimization basin $(select_idx)/$(n_catchments)")
        recorder = []
        last_time = time()
        callback_func!(state, l) = begin
            current_time = time()
            delta_time = current_time - last_time
            push!(recorder, (iter=state.iter, loss=l, time=now(), delta_time=delta_time, params=state.u))
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
        loss_history = recorder_df.loss
        cumulative_time = cumsum(recorder_df.delta_time)
        @save result_file params loss_history train_kge test_kge train_loss test_loss elapsed_seconds cumulative_time recorder_df

        push!(metrics_df, (basin=select_idx, train_kge=train_kge, test_kge=test_kge, train_loss=train_loss, test_loss=test_loss))
        all_recorders[select_idx] = recorder
        all_params[select_idx] = params
        @info "basin metrics" basin = select_idx train_kge = train_kge test_kge = test_kge train_loss = train_loss test_loss = test_loss
    end

    GC.gc()
end

@info "batch metrics" metrics_df
# save aggregate metrics and parameters across basins
batch_metrics_file = joinpath(result_dir, "batch_metrics.jld2")
@save batch_metrics_file metrics_df all_params all_recorders