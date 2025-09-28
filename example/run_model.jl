# Import required packages
using CSV, DataFrames, JLD2
using ComponentArrays
using Distributions
using HydroModels, HydroModelTools
using DataInterpolations
using Symbolics: tosymbol
using Accessors
using Plots
include("../src/HydroModelLibrary.jl")

# Data Loading: Read input data from CSV and prepare time index
df = CSV.read("data/marrmot/3604000.csv", DataFrame)
tidx = rem.(0:length(df[!, "prec"])-1, 365) .+ 1
input = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"], Tidx=tidx)

test_input = stack(input, dims=1)

# Model Setup: Load and initialize the Penman model
model_nm = "sacramento"
model_module = HydroModelLibrary.load_model(Symbol(model_nm), reload=true)
model = model_module.model
model_variables = model_module.model_variables
model_parameters = model_module.model_parameters

# Parameter Preparation: Define parameter names and bounds, initialize parameters
init_params = HydroModelLibrary.get_random_params(model_nm)
init_params = JLD2.load("cache/calibrate/$model_nm/sol.jld2", "opt_params")
@info "init_params: $init_params"

# Model Input/Output Setup: Define model inputs, states, and outputs
model_input_names = HydroModels.get_input_names(model)
model_state_names = HydroModels.get_state_names(model)
model_output_names = HydroModels.get_output_names(model)
init_states = NamedTuple{Tuple(model_state_names)}(zeros(length(model_state_names))) |> ComponentVector
@info "Input variables: $(model_input_names)"

# Model Execution: Prepare input array and run the model
input_arr = stack(input[model_input_names], dims=1)
config = (; solver=HydroModelTools.ODESolver(), interp=LinearInterpolation)
result = model(input_arr, init_params, initstates=init_states, config=config)

# Output Processing: Convert results to DataFrame and compute summaries
model_output_names = vcat(model_state_names, model_output_names) |> Tuple
output_df = DataFrame(NamedTuple{model_output_names}(eachslice(result, dims=1)))
r2_func(y, y_hat) = sum((y .- y_hat) .^ 2) ./ sum((y .- mean(y)) .^ 2)
@info "sum of flow" df[!, "flow"] |> sum
@info "sum of Qt" output_df[!, "Qt"] |> sum
@info "sum of prec" df[!, "prec"] |> sum
@info 1 - r2_func(df[!, "flow"], output_df[!, "Qt"])

# Visualization: Plot observed flow and model output
plot(df[750:1000, "flow"], label="flow")
plot!(output_df[750:1000, "Qt"], label="Qt")