using CSV, DataFrames
include("../src/HydroModelLibrary.jl")

df = CSV.read("data/marrmot/3604000.csv", DataFrame)
tidx = rem.(0:length(df[!, "prec"])-1, 365) .+ 1
input = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"], Tidx=tidx)
flow_vec = df[!, "flow"]
record = Dict{Symbol, Any}()
r2_func(y, y_hat) = sum((y .- y_hat) .^ 2) ./ sum((y .- mean(y)) .^ 2)
for model_nm in HydroModelLibrary.AVAILABLE_MODELS
    load_path = "cache/calibrate/$model_nm"
    output_path = "$load_path/output.csv"
    tmp_output_df = CSV.read(output_path, DataFrame)
    r2 = 1 - r2_func(flow_vec, tmp_output_df[!, :Qt])
    record[model_nm] = r2
end
# 打印r2小于0.5的模型
for model_nm in keys(record)
    if record[model_nm] < 0.5
        @info model_nm, record[model_nm]
    end
end