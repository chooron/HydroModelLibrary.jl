# HydroModelLibrary.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://chooron.github.io/HydroModelLibrary.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

`HydroModelLibrary.jl` is a hydrological model library built on top of [HydroModels.jl](https://github.com/chooron/HydroModels.jl). It collects multiple conceptual rainfall-runoff models, mainly from [MARRMoT](https://github.com/wknoben/MARRMoT), and provides example workflows for simulation, calibration, and result inspection.

## Features

- Multiple hydrological model implementations that can be loaded directly
- Unified access to model parameters, parameter bounds, states, inputs, and outputs
- Example scripts under `example/` for running models, batch calibration, and checking saved results

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

### Run a single model

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

### Calibrate `hbv_edu` on CAMELS PET v2

The main batch-calibration example is [`example/calibrate_hbv_camels_petv2.jl`](example/calibrate_hbv_camels_petv2.jl). It calibrates `hbv_edu` basin by basin on the CAMELS PET v2 dataset and writes per-basin and batch summary results to disk.

Current script behavior:

- Working directory should be `example/`
- Data path is read from `.env`
- Environment variable name is `CAMESL_DATASET_PATH`
- Objective function defaults to `inverse KGE`
- Training period is `1989-01-01` to `1998-12-31`
- Test period is `1999-01-01` to `2009-12-31`
- Warm-up length is `365` days
- Output directory is `example/calibration_results/hbv_invkge/`

Configure the dataset path in `example/.env`:

```env
CAMESL_DATASET_PATH=G:\Dataset\camels_dataset_petv2.npz
```

Run from the `example/` directory:

```powershell
cd E:\JlCode\HydroModelLibrary\example
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. calibrate_hbv_camels_petv2.jl
```

## Input Data

`example/calibrate_hbv_camels_petv2.jl` expects an `.npz` file with at least:

- `forcing`, shaped as `(time, basin, variable)`
- `target`, shaped as `(time, basin)`

The script interprets `forcing[:, :, 1:3]` as:

- `P`: precipitation
- `T`: temperature
- `Ep`: potential evapotranspiration

It also builds `Tidx = dayofyear.(dates)` as an additional input required by `hbv_edu`.

## Output Files

For each processed basin, the script writes:

```text
example/calibration_results/hbv_invkge/basin_<idx>.jld2
```

Each per-basin result file stores:

- `params`
- `loss_history`
- `train_kge`
- `test_kge`
- `train_loss`
- `test_loss`
- `elapsed_seconds`
- `cumulative_time`
- `recorder_df`

After the loop finishes, the script also writes:

```text
example/calibration_results/hbv_invkge/batch_metrics.jld2
```

This batch file stores:

- `metrics_df`
- `all_params`
- `all_recorders`

## Notes

- The current loop is `for select_idx in 2:n_catchments`, so `basin_1` is skipped by default
- If a basin result file already exists, the script skips recalibration for that basin and reuses the saved outputs
- `forcing` and `target` must have the same time length
- Both train and test splits must be longer than the `365`-day warm-up

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
