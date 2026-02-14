using JLD2
using DataFrames
using Statistics
using ComponentArrays

models = Dict(
	"collie1" => "example/calibration_results/collie1_kge",
	"xaj" => "example/calibration_results/xaj_kge",
	"hbv" => "example/calibration_results/hbv_kge",
	"hymod" => "example/calibration_results/hymod_kge",
)

function safe_median(values)
	finite_values = filter(isfinite, values)
	return isempty(finite_values) ? NaN : median(finite_values)
end

function model_stats(model_dir::AbstractString)
	files = filter(f -> startswith(f, "basin_") && endswith(f, ".jld2"), readdir(model_dir))
	iter_times = Float64[]
	min_updates = Float64[]
	for file in files
		info = JLD2.load(joinpath(model_dir, file))
		iter_time = info["elapsed_seconds"] / length(info["cumulative_time"])
		push!(iter_times, iter_time)
		push!(min_updates, 60 / iter_time)
	end
	return safe_median(iter_times), safe_median(min_updates), length(files)
end

rows = DataFrame(model = String[], median_iter_time = Float64[], median_min_update = Float64[], basin_count = Int[])
for (name, dir) in models
	median_iter, median_min_update, n = model_stats(dir)
	push!(rows, (name, median_iter, median_min_update, n))
end

rows