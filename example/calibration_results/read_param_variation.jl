using JLD2
using ComponentArrays
using DataFrames
using Plots

basin_id = 100
model_name = "xaj"
result = JLD2.load(joinpath("example", "calibration_results", "$(model_name)_kge", "basin_$(basin_id).jld2"))
record_df = result["recorder_df"]
record_params_list = record_df[!, "params"]
param_arr = stack(record_params_list, dims=1)

