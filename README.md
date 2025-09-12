# HydroModelLibrary.jl

## Overview
HydroModelLibrary.jl is a hydrological modeling library built upon [HydroModels.jl](https://github.com/chooron/HydroModels.jl). It incorporates models primarily sourced from [MARRMoT](https://github.com/wknoben/MARRMoT) and includes additional models such as ExpHydro. All models provided in HydroModelLibrary.jl have been validated and tested for direct use.

## Usage
To load and use a model, you can utilize the following code snippet:

```julia
model_module = HydroModelLibrary.load_model(:penman, reload=true)
model = model_module.model
```

## Features
- **Model Implementation**: The library includes a variety of hydrological models, with ongoing efforts to implement runoff computation modules in future updates.
- **Examples**: The project provides example scripts located at:
  - `example/run_model.jl`: Demonstrates model computation.
  - `example/optimal_hydromodels.jl`: Implements parameter optimization methods for hydrological models.
  
Users can copy these scripts and adapt them with their own data for calibration purposes.

## Notes
- [HydroModels.jl](https://github.com/chooron/HydroModels.jl) and its ecosystem are under active development. Ensure all dependencies are updated to the latest versions.
- If you encounter issues, please submit them via the project's [issue tracker](https://github.com/chooron/HydroModels.jl/issues).