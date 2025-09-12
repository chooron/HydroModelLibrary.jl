using CSV, DataFrames
using ComponentArrays
using Distributions
using HydroModels, HydroModelTools
using DataInterpolations
using Symbolics: tosymbol
using Accessors
using Plots
include("../src/HydroModelLibrary.jl")

df = CSV.read("data/marrmot/3604000.csv", DataFrame)
tidx = rem.(0:length(df[!, "prec"])-1, 365) .+ 1

input = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"], Tidx=tidx)
model_module = HydroModelLibrary.load_model(:penman, reload=true)
model = model_module.model
model_variables = model_module.model_variables
model_parameters = model_module.model_parameters
model_params_names = tosymbol.(model_parameters)
model_params_bounds = NamedTuple{Tuple(model_params_names)}(HydroModels.getbounds.(model_parameters))
random_param_values = map(zip(model_params_names, model_params_bounds)) do (param_name, param_bound)
    param_name => rand(Uniform(param_bound[1], param_bound[2]))
end |> NamedTuple
init_params = ComponentVector(params=random_param_values)
init_params = JLD2.load("cache/calibrate/penman/sol.jld2", "opt_params")
@info "init_params: $init_params"
model_input_names = HydroModels.get_input_names(model)
model_state_names = HydroModels.get_state_names(model)
model_output_names = HydroModels.get_output_names(model)
init_states = NamedTuple{Tuple(model_state_names)}(zeros(length(model_state_names))) |> ComponentVector
@info "Input variables: $(model_input_names)"
input_arr = stack(input[model_input_names], dims=1)
config = (; solver=HydroModelTools.ODESolver(), interp=LinearInterpolation)
result = model(input_arr, init_params, initstates=init_states, config=config)
model_output_names = vcat(model_state_names, model_output_names) |> Tuple
output_df = DataFrame(NamedTuple{model_output_names}(eachslice(result, dims=1)))
@info "sum of flow" df[!, "flow"] |> sum
@info "sum of Qt" output_df[!, "Qt"] |> sum
@info "sum of prec" df[!, "prec"] |> sum
@info r2_func(df[!, "flow"], output_df[!, "Qt"])
plot(df[!, "flow"], label="flow")
plot!(output_df[!, "Qt"], label="Qt")

plot(df[500:900, "flow"], label="flow")
plot!(output_df[500:900, "Qt"], label="Qt")
# plot!(output_df[500:700, "Uq"], label="q")