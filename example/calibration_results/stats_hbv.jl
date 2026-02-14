using JLD2
using DataFrames
using Statistics

const BASE_DIR = "example/calibration_results/hbv_invkge"

function collect_hbv_kge(basin_range=2:559)
    df = DataFrame(
        basin=Int[],
        train_kge=Vector{Union{Missing,Float64}}(),
        test_kge=Vector{Union{Missing,Float64}}(),
    )

    for i in basin_range
        path = joinpath(BASE_DIR, "basin_$(i).jld2")
        isfile(path) || continue
        data = JLD2.load(path)
        train_val = get(data, "train_kge", get(data, :train_kge, missing))
        test_val = get(data, "test_kge", get(data, :test_kge, missing))
        push!(df, (basin=i, train_kge=train_val, test_kge=test_val))
    end

    return df
end


df = collect_hbv_kge()

train_vals = filter(isfinite, collect(skipmissing(df.train_kge)))
test_vals = filter(isfinite, collect(skipmissing(df.test_kge)))

diff_count = sum((ismissing(tr) || ismissing(ts) || !isfinite(tr) || !isfinite(ts)) ? 0 : (tr - ts > 0.2) for (tr, ts) in zip(df.train_kge, df.test_kge))

println(df)

if !isempty(train_vals)
    println("train_kge mean=$(mean(train_vals)) median=$(median(train_vals))")
end

if !isempty(test_vals)
    println("test_kge mean=$(mean(test_vals)) median=$(median(test_vals))")
end

println("count(train_kge - test_kge > 0.5) = $(diff_count)")

"""
train_kge mean=0.8305628205192583 median=0.8548066044989484
test_kge mean=0.6468509462583718 median=0.7296854401459514
count(train_kge - test_kge > 0.5) = 163
"""