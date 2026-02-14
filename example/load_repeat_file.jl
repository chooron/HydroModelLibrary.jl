using JLD2
using DataFrames
using ComponentArrays
using NPZ

const BASIN_IDX = 7
const REPEAT_COUNT = 10

model_name = "hbv96"
result_dir = joinpath(@__DIR__, "calibration_results", "$(model_name)_repeat", "basin_$(BASIN_IDX)")
mkpath(result_dir)
params_arr_list = Vector{Array{Float64,2}}()

for run_idx in 1:REPEAT_COUNT
    result_file = joinpath(result_dir, "basin_$(BASIN_IDX)_run_$(run_idx).jld2")
    tmp_result = JLD2.load(result_file)
    tmp_record_df = tmp_result["recorder_df"]
    params_list = tmp_record_df[:, "params"]
    params_arr = stack(params_list, dims=1)  # iterations × parameters
    push!(params_arr_list, params_arr)
end

sizes = map(size, params_arr_list)
@assert length(unique(sizes)) == 1 "Parameter trajectories have mismatched shapes across runs"

# 3D tensor: runs × iterations × parameters
param_tensor = cat(params_arr_list...; dims=3) |> x -> permutedims(x, (3, 1, 2))
param_tensor[:,:,2]
output_npy = joinpath(result_dir, "$(model_name)_runs_iters_params.npy")
NPZ.npzwrite(output_npy, param_tensor)