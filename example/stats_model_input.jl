include("../src/HydroModelLibrary.jl")

model_input_dict = Dict{Symbol, Any}()

for model_nm in HydroModelLibrary.AVAILABLE_MODELS
    model_module = HydroModelLibrary.load_model(Symbol(model_nm))
    model = model_module.model
    model_input_names = HydroModels.get_input_names(model)
    model_input_dict[model_nm] = model_input_names
end

for model_nm in keys(model_input_dict)
    if length(model_input_dict[model_nm]) == 3
        @info model_nm
    end
end