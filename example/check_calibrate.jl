using JLD2
using ComponentArrays
using Plots

results1 = JLD2.load(joinpath(@__DIR__, "calibration_results", "xaj_kge", "basin_1.jld2"))
results2 = JLD2.load(joinpath(@__DIR__, "calibration_results", "collie1_kge", "basin_1.jld2"))
plot(results1["cumulative_time"], results1["loss_history"])
plot!(results2["cumulative_time"], results2["loss_history"])