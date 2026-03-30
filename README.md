# HydroModelLibrary.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://chooron.github.io/HydroModelLibrary.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

`HydroModelLibrary.jl` is a hydrological model library built on top of [HydroModels.jl](https://github.com/chooron/HydroModels.jl). It collects multiple conceptual rainfall-runoff models, mainly from [MARRMoT](https://github.com/wknoben/MARRMoT), and provides example workflows for simulation, calibration, and result inspection.

## Features

- Multiple hydrological model implementations that can be loaded directly
- Unified access to model parameters, parameter bounds, states, inputs, and outputs
- Example scripts under `example/` for model runs, batch calibration, and result inspection

## Installation

Install the registered package:

```julia
using Pkg
Pkg.add("HydroModelLibrary")
```

Or install the latest development version from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/chooron/HydroModelLibrary.jl")
```

## Quick Start

### Run a Single Model

```julia
using CSV
using DataFrames
using ComponentArrays
using DataInterpolations
using HydroModels
using HydroModelLibrary

df = CSV.read("data/marrmot/3604000.csv", DataFrame)
input_data = (P=df.prec, Ep=df.pet, T=df.temp)

model_module = HydroModelLibrary.load_model(:sacramento)
model = model_module.model
model_input_names = HydroModels.get_input_names(model)
model_state_names = HydroModels.get_state_names(model)

init_params = HydroModelLibrary.get_random_params("sacramento")
init_states = NamedTuple{Tuple(model_state_names)}(zeros(length(model_state_names))) |> ComponentVector
input_matrix = stack(input_data[model_input_names], dims=1)

config = (; solver=HydroModels.ODESolver, interp=LinearInterpolation)
results = model(input_matrix, init_params; initstates=init_states, config=config)

println(size(results))
```

### `calibrate_hbv_camels_petv2.jl` Walkthrough

The following sections explain [`example/calibrate_hbv_camels_petv2.jl`](example/calibrate_hbv_camels_petv2.jl) in execution order. The code blocks below are the script content, split into logical sections and explained step by step.

#### 1. Activate the Environment and Import Dependencies

This section activates the Julia environment under `example/`, loads the required packages, reads configuration from `.env`, and includes the local package source.

```julia
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
```

This script assumes:

- The working directory is `example/`
- `.env` defines `CAMESL_DATASET_PATH`

Example `.env` entry:

```env
CAMESL_DATASET_PATH=G:\Dataset\camels_dataset_petv2.npz
```

Run it with:

```powershell
cd E:\JlCode\HydroModelLibrary\example
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. calibrate_hbv_camels_petv2.jl
```

#### 2. Load the CAMELS PET v2 Dataset

This section reads the `.npz` dataset and determines how many basins are available.

```julia
# load data
data = npzread(ENV["CAMESL_DATASET_PATH"])
n_catchments = size(data["forcing"], 2)
```

The input file must contain at least:

- `forcing`, shaped as `(time, basin, variable)`
- `target`, shaped as `(time, basin)`

Later in the script, `forcing[:, :, 1:3]` is interpreted as `P`, `T`, and `Ep`.

#### 3. Build the Date Index and Define Train/Test Splits

This section creates the full date range, trims it to the available data length, defines training and testing periods, and sets the warm-up length.

```julia
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
```

This part enforces several assumptions:

- `forcing` and `target` must have the same number of time steps
- The training period is fixed to `1989-01-01` through `1998-12-31`
- The test period is fixed to `1999-01-01` through `2009-12-31`
- Both splits must be longer than the `365`-day warm-up

#### 4. Define the Calibration Configuration

This section sets the model name, optimization budget, and whether to use inverse KGE.

```julia
max_iter = 10000
model_nm = "hbv_edu"
const USE_INVERSE_KGE = true
const INVKGE_EPS_FACTOR = 0.01
const INVKGE_MIN_EPS = 1e-3
const INVKGE_STABILITY_EPS = 1e-5
```

The script is currently hard-coded to calibrate `hbv_edu`. It also enables `inverse KGE` by default, which emphasizes low-flow behavior by transforming discharge before computing KGE.

#### 5. Load the Model and Extract Metadata

This section loads `hbv_edu` from `HydroModelLibrary` and extracts parameter names, input names, state names, output names, and parameter bounds.

```julia
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
```

This metadata is built once and then reused for all basins.

#### 6. Prepare Shared Bounds, Solver Settings, and Output Containers

This section converts parameter bounds into vectors for the optimizer, defines the solver configuration, creates the output directory, and initializes result containers.

```julia
# bounds and config shared by all basins
lb_list = first.(model_params_bounds |> collect)
ub_list = last.(model_params_bounds |> collect)
config = (solver=HydroModels.ODESolver, interp=LinearInterpolation)
result_dir = joinpath(@__DIR__, "calibration_results", "hbv_invkge")
mkpath(result_dir)

metrics_df = DataFrame(basin=Int[], train_kge=Float64[], test_kge=Float64[], train_loss=Float64[], test_loss=Float64[])
all_recorders = Vector{Any}(undef, n_catchments)
all_params = Vector{Any}(undef, n_catchments)
```

Outputs are written under:

```text
example/calibration_results/hbv_invkge/
```

The three main containers are:

- `metrics_df` for per-basin train/test metrics
- `all_recorders` for per-basin iteration logs
- `all_params` for per-basin calibrated parameters

#### 7. Define KGE and Inverse KGE

This section implements both standard KGE and inverse KGE. Both functions return `(loss, kge)`.

```julia
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
```

Two implementation details matter here:

- Evaluation always skips the first `365` warm-up days
- `sim[end, eval_idxs]` assumes the last model output row is streamflow

#### 8. Run Basin-by-Basin Calibration

This is the main loop. For each basin, the script checks whether a result file already exists. If it does, the script loads the saved metrics. Otherwise, it prepares basin-specific data, samples a random initial parameter vector, runs optimization, evaluates train/test performance, and saves the result.

```julia
for select_idx in 2:n_catchments
    result_file = joinpath(result_dir, "basin_$(select_idx).jld2")

    # If the basin result file already exists, load saved metrics and skip recalibration
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
        @info "basin already processed; skip" basin=select_idx
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
```

This loop can be understood in four steps:

- Check whether `basin_<idx>.jld2` already exists and skip recalibration if it does
- Build basin-specific forcing and target vectors, then assemble the `(P, T, Ep, Tidx)` input tuple
- Draw a random initial parameter set within bounds and optimize the training loss with `BBO_adaptive_de_rand_1_bin_radiuslimited()`
- Save parameters, loss history, train metrics, test metrics, and the iteration recorder

There are also three concrete details worth calling out:

- The loop is `for select_idx in 2:n_catchments`, so `basin_1` is skipped by default
- The objective function uses only the training split loss for optimization
- Even though the optimizer is black-box, the script still wraps the objective with `Optimization.AutoForwardDiff()`

#### 9. Save the Batch Summary

After the basin loop finishes, the script writes a final summary file.

```julia
@info "batch metrics" metrics_df
# save aggregate metrics and parameters across basins
batch_metrics_file = joinpath(result_dir, "batch_metrics.jld2")
@save batch_metrics_file metrics_df all_params all_recorders
```

The final outputs are:

- Per-basin files: `example/calibration_results/hbv_invkge/basin_<idx>.jld2`
- Batch summary: `example/calibration_results/hbv_invkge/batch_metrics.jld2`

#### 10. Practical Notes About This Script

- The environment variable name is literally `CAMESL_DATASET_PATH`, so `.env` and documentation must use that spelling
- The training and testing periods are fixed in the script, not passed in externally
- Evaluation always excludes the first `365` days as warm-up
- Existing result files are automatically reused, which makes interrupted runs resumable
- The current implementation starts from basin 2, not basin 1

## Example Scripts

Other scripts under `example/` include:

- `run_model.jl`
- `load_sample_data_example.jl`
- `calibrate_hymod_camels_petv2.jl`
- `calibrate_xaj_camels_petv2.jl`
- `calibrate_collie1_camels_petv2.jl`
- `check_calibrate.jl`

## Ecosystem

This package relies mainly on:

- [HydroModels.jl](https://github.com/chooron/HydroModels.jl)
- [SciML](https://sciml.ai/)
- [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl)
- [Optimization.jl](https://github.com/SciML/Optimization.jl)
- [JLD2.jl](https://github.com/JuliaIO/JLD2.jl)

## Citation

If you use `HydroModelLibrary.jl` in research, please cite:

```bibtex
@software{hydromodellibrary_jl,
  author = {Jing Xin},
  title = {HydroModelLibrary.jl: A Library of Hydrological Models},
  year = {2024},
  url = {https://github.com/chooron/HydroModelLibrary.jl}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Author

**Jing Xin** ([@chooron](https://github.com/chooron))
