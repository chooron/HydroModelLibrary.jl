module collie1
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables S [description = "current storage in the soil moisture", unit = "mm"]
@variables P [description = " precipitation input", unit = "mm/d"]
@variables Ea [description = "estimated based on the current storage S", unit = "mm/d"]
@variables Smax [description = "the maximum soil moisture storage", unit = "mm"] 
@variables Ep [description = "potential evapotranspiration", unit = "mm/d"]
@variables Qt [description = "saturation excess overland flow", unit = "mm/d"]
model_variables = [S, P, Ea, Smax, Ep, Qt]

# Model parameters
@parameters Smax [description = "maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
model_parameters = [Smax]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ea ~  min(clamp(S /  Smax, 0, 1) * Ep, S)
        @hydroflux Qt ~ step_func(S - Smax) * P
    end
    dfluxes = begin
        @stateflux S ~ P  - Ea - Qt
    end
end

model = @hydromodel :collie1 begin
    bucket1
end

end
