using CSV
using JLD2
using DataFrames
using Statistics

models = Dict(
	"collie1" => "example/calibration_results/collie1_invkge",
	"xaj" => "example/calibration_results/xaj_invkge",
	"hbv" => "example/calibration_results/hbv_invkge",
	"hymod" => "example/calibration_results/hymod_invkge",
)

function epoch_medians(model_dir::AbstractString)
	files = filter(f -> startswith(f, "basin_") && endswith(f, ".jld2"), readdir(model_dir))
	time_lists = Dict{Int, Vector{Float64}}()
	loss_lists = Dict{Int, Vector{Float64}}()

	for file in files
		info = JLD2.load(joinpath(model_dir, file))
		cumulative_time = get(info, "cumulative_time", nothing)
		loss_history = get(info, "loss_history", nothing)
		cumulative_time === nothing && continue
		loss_history === nothing && continue

		n = min(length(cumulative_time), length(loss_history))
		for epoch in 1:n
			ct = cumulative_time[epoch]
			ls = loss_history[epoch]
			if isfinite(ct) && isfinite(ls)
				push!(get!(time_lists, epoch, Float64[]), ct)
				push!(get!(loss_lists, epoch, Float64[]), ls)
			end
		end
	end

	epochs = sort!(collect(keys(time_lists)))
	df = DataFrame(epoch = Int[], median_time = Float64[], median_loss = Float64[])
	for epoch in epochs
		tlist = filter(isfinite, time_lists[epoch])
		llist = filter(isfinite, loss_lists[epoch])
		(isempty(tlist) || isempty(llist)) && continue
		push!(df, (epoch, median(tlist), median(llist)))
	end

	return df, length(files)
end

results = Dict{String, DataFrame}()
for (name, dir) in models
	df, basin_count = epoch_medians(dir)
	results[name] = df
	out_path = joinpath(dir, "$(name)_invkge_loss_curve_medians.csv")
	CSV.write(out_path, df)
	@info "saved loss curve stats" model = name basins = basin_count epochs = nrow(df) path = out_path
end

results