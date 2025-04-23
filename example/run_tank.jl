using CSV, DataFrames
using ComponentArrays
using HydroModels
using HydroModelTools
using Distributions
using ModelingToolkit
using DataInterpolations
using Plots
include("../src/HydroModelLibrary.jl")

df = CSV.read("data/marrmot/3604000.csv", DataFrame)

input = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"])
model = HydroModelLibrary.tank.model
param_bounds = getbounds.(get_params(model))
random_param_values = map(param_bounds) do param_bound
    rand(Uniform(param_bound[1], param_bound[2]))
end
init_params = ComponentVector(params=NamedTuple{Tuple(get_param_names(model))}(random_param_values))

t2 = init_params.params.f2 * init_params.params.st         # Threshold before surface runoff 2 starts [mm]
t1 = t2 + init_params.params.f1 * (init_params.params.st - t2) # Threshold before surface runoff 1 starts [mm]
t3 = init_params.params.f3 * (init_params.params.st - t1)    # Threshold before intermediate runoff starts [mm]
t4 = init_params.params.st - t1 - t3      # Threshold before sub-base runoff starts [mm]

init_states = NamedTuple{Tuple(get_state_names(model))}(ones(length(get_state_names(model)))) |> ComponentVector
@info "Input variables: $(HydroModels.get_input_names(model))"
input_arr = stack(input[HydroModels.get_input_names(model)], dims=1)
config = (; solver=HydroModelTools.ODESolver(), interp=LinearInterpolation)
result = model(input_arr, init_params, initstates=init_states, config=config)
model_output_names = vcat(get_state_names(model), get_output_names(model)) |> Tuple
output_df = DataFrame(NamedTuple{model_output_names}(eachslice(result, dims=1)))

plot(output_df[!, "Qt"], label="Q_hat")
plot!(df[!, "flow"], label="Q")

plot(output_df[!, "S1"], label="S1")
plot!(output_df[!, "S2"], label="S2")
plot!(output_df[!, "S3"], label="S3")
plot!(output_df[!, "S4"], label="S4")

