# HydroModelLibrary.jl

[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://chooron.github.io/HydroModelLibrary.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## Overview

A hydrological modeling library built upon [HydroModels.jl](https://github.com/chooron/HydroModels.jl), incorporating models primarily sourced from [MARRMoT](https://github.com/wknoben/MARRMoT) and additional models like ExpHydro. All models are validated and ready for direct use.

## ✨ Key Features

- **Model Implementation**: Diverse hydrological models, with ongoing efforts for runoff computation modules.
- **Examples**: Comprehensive examples for model computation and parameter optimization.

## 🚀 Installation

```julia
using Pkg
Pkg.add("HydroModelLibrary")
```

Or install from GitHub for the latest development version:

```julia
Pkg.add(url="https://github.com/chooron/HydroModelLibrary.jl")
```

## 📚 Quick Start

### Example 1: Running a Hydrological Model

Demonstrates how to load a model, prepare input data, and run simulations.

```julia
using HydroModelLibrary
using CSV, DataFrames, ComponentArrays
using HydroModels, HydroModelTools
using DataInterpolations # For input interpolation

# Load sample data (assuming 'data/marrmot/3604000.csv' exists)
df = CSV.read("data/marrmot/3604000.csv", DataFrame)
input_data = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"])

# Load a specific model, e.g., "sacramento"
model_name = "sacramento"
model_module = HydroModelLibrary.load_model(Symbol(model_name), reload=true)
model = model_module.model
model_input_names = HydroModels.get_input_names(model)
model_state_names = HydroModels.get_state_names(model)

# Prepare initial parameters and states
init_params = HydroModelLibrary.get_random_params(model_name)
init_states = NamedTuple{Tuple(model_state_names)}(zeros(length(model_state_names))) |> ComponentVector

# Prepare input array
input_arr = stack(input_data[model_input_names], dims=1)

# Run the model
config = (; solver=HydroModelTools.ODESolver(), interp=LinearInterpolation)
results = model(input_arr, init_params, initstates=init_states, config=config)

println("Model simulation completed. Results dimensions: ", size(results))
# Further processing and visualization of results can be done here.
```

### Example 2: Gradient-Based Parameter Optimization

Illustrates how to set up and run parameter optimization for a hydrological model.

```julia
using HydroModelLibrary
using CSV, DataFrames, ComponentArrays
using HydroModels, HydroModelTools
using DataInterpolations
using Optimization, OptimizationBBO # For optimization algorithms
using Distributions
using Symbolics # To work with symbolic parameters

# Load sample data and observed flow
df = CSV.read("data/marrmot/3604000.csv", DataFrame)
input_data = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"])
observed_flow = df[!, "flow"]

# Load a specific model, e.g., "penman"
model_module = HydroModelLibrary.load_model(:penman)
model = model_module.model
model_parameters = model_module.model_parameters
model_params_names = tosymbol.(model_parameters)
model_input_names = HydroModels.get_input_names(model)

# Define objective function (e.g., R2 loss)
function obj_func(p, _)
    params = ComponentVector(p, getaxes(HydroModelLibrary.get_random_params("penman"))) # Reconstruct ComponentVector
    input_matrix = stack(input_data[model_input_names], dims=1)
    config = (solver=HydroModelTools.ODESolver(), interp=LinearInterpolation)
    predictions = model(input_matrix, params; config=config)
    predicted_flow = predictions[end, :] # Assuming flow is the last output

    # Calculate R2 loss (simplified for demonstration)
    y_mean = mean(observed_flow[365:end]) # Warm-up period
    sum_sq_res = sum((predicted_flow[365:end] .- observed_flow[365:end]) .^ 2)
    sum_sq_tot = sum((observed_flow[365:end] .- y_mean) .^ 2)
    return sum_sq_res / sum_sq_tot # 1 - R2 for minimization
end

# Prepare initial parameters and bounds
init_params = HydroModelLibrary.get_random_params(:penman)
model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModels.getbounds.(model_parameters))
lb = first.(model_params_bounds |> collect)
ub = last.(model_params_bounds |> collect)

# Setup and run optimization
optf = Optimization.OptimizationFunction(obj_func, Optimization.AutoForwardDiff())
optprob = Optimization.OptimizationProblem(optf, Vector(init_params), lb=lb, ub=ub)
sol = Optimization.solve(optprob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=1000)

println("Optimization completed. Optimal parameters: ", ComponentVector(sol.u, getaxes(init_params)))
```

## 🔗 Ecosystem Integration

HydroModelLibrary.jl integrates seamlessly with the Julia ecosystem, leveraging:
- **[HydroModels.jl](https://github.com/chooron/HydroModels.jl)**: The foundational modeling framework.
- **[SciML](https://sciml.ai/)**: For ODE solving, sensitivity analysis, and optimization.
- **[ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl)**: For named parameter arrays.
- **[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)**: For symbolic computation.
- **[Optimization.jl](https://github.com/SciML/Optimization.jl)**: For advanced optimization algorithms.

## 📝 Citation

If you use HydroModelLibrary.jl in your research, please cite:

```bibtex
@software{hydromodellibrary_jl,
  author = {Jing Xin},
  title = {HydroModelLibrary.jl: A Library of Hydrological Models},
  year = {2024},
  url = {https://github.com/chooron/HydroModelLibrary.jl}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author

**Jing Xin** ([@chooron](https://github.com/chooron))